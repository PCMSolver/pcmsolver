#include "PWCSolver.hpp"

#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>
//#include "Getkw.h"

extern "C"
{
#include "vector3.h"
#include "sparse2.h"
#include "intvector.h"
#include "basis.h"
#include "WEM.h"
#include "read_points.h"
#include "vector2.h"
#include "interpolate.h"
#include "topology.h"
#include "kern.h"
#include "compression.h"
#include "postproc.h"
#include "WEMRHS.h"
#include "WEMPCG.h"
#include "WEMPGMRES.h"
#include "dwt.h"
#include "cubature.h"
#include "gauss_square.h"
#include "constants.h"
#include "energy.h"
}

#include "Cavity.hpp"
#include "GreensFunction.hpp"
#include "WaveletCavity.hpp"

static GreensFunction * gf;

static double SLInt(vector3 x, vector3 y)
{  
	return(1.0l/sqrt((x.x-y.x)*(x.x-y.x)+(x.y-y.y)*(x.y-y.y)+(x.z-y.z)*(x.z-y.z)));
}


static double SLExt(vector3 x, vector3 y)
{
	double 		r = sqrt((x.x-y.x)*(x.x-y.x)+(x.y-y.y)*(x.y-y.y)+(x.z-y.z)*(x.z-y.z));
	return 1.0l/(r*78.39l);
}


static double DLUni(vector3 x, vector3 y, vector3 n_y)
{  
	vector3		c;
	double		r;
	Eigen::Vector3d grad, dir;
	c.x = x.x-y.x;
	c.y = x.y-y.y;
	c.z = x.z-y.z;
	r = sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
	grad(0) = c.x/(r*r*r);
	grad(1) = c.y/(r*r*r);
	grad(2) = c.z/(r*r*r);
	dir << n_y.x, n_y.y, n_y.z;
	return (c.x*n_y.x+c.y*n_y.y+c.z*n_y.z)/(r*r*r);
}

static double SingleLayer (vector3 x, vector3 y) 
{
	Eigen::Vector3d vx(x.x, x.y, x.z);
	Eigen::Vector3d vy(y.x, y.y, y.z);
	Eigen::Vector3d foo = Eigen::Vector3d::Zero();
	double value = gf->evaluate(foo, vx, foo, vy)(0);
	return value;
}

static double DoubleLayer (vector3 x, vector3 y, vector3 n_y) 
{
	Eigen::Vector3d vx(x.x, x.y, x.z);
	Eigen::Vector3d vy(y.x, y.y, y.z);
	Eigen::Vector3d vn_y(n_y.x, n_y.y, n_y.z);
	Eigen::Vector3d foo = Eigen::Vector3d::Zero();
	double value = gf->evaluate(foo, vx, vn_y, vy)(1);
	return value;
}

void PWCSolver::initPointers()
{
	elementTree = NULL;
	waveletList = NULL;
}


/*PWCSolver::PWCSolver(const Section & solver) : WEMSolver(solver) {
	initPointers();
	setSolverType("Wavelet");
}*/

PWCSolver::~PWCSolver(){
	if(elementTree != NULL) free_elementlist(&elementTree, nPatches,nLevels);
	if(waveletList != NULL) free_waveletlist(&waveletList, nPatches,nLevels);
	if(T_ != NULL)          free_interpolate(&T_,nPatches,nLevels);
}

void PWCSolver::initInterpolation() {
	init_interpolate(&T_, pointList, nPatches, nLevels);
	nNodes = gennet(&nodeList, &elementList, pointList, 
					nPatches, nLevels);
}

void PWCSolver::constructWavelets(){
	generate_elementlist(&elementTree, nodeList, elementList, nPatches, nLevels);
	generate_waveletlist(&waveletList, elementTree, nPatches, nLevels);
	set_quadrature_level(waveletList, elementTree, nPatches, nLevels);
	simplify_waveletlist(waveletList, elementTree, nPatches, nLevels);
	complete_elementlist(waveletList, elementTree, nPatches, nLevels);
}	

void PWCSolver::constructSi() {
	double factor = 0;
	double epsilon = 0;
	switch (integralEquation) {
	case FirstKind:
		epsilon = greenOutside->getDielectricConstant();
		factor = - 2 * M_PI * (epsilon + 1) / (epsilon - 1);
		break;
	case SecondKind:
		throw std::runtime_error("Second Kind not yet implemented"); //careful to the double layer sign when implementing it....
		break;
	case Full:
		factor = 2 * M_PI;
		break;
	default:
		exit(-1);
	}
	gf = greenInside;
	apriori1_ = compression(&S_i_, waveletList, elementTree, nPatches, nLevels);
	WEM(&S_i_, waveletList, elementTree, T_, nPatches, nLevels, 
			SingleLayer, DoubleLayer, factor);
	aposteriori1_ = postproc(&S_i_, waveletList, elementTree, nPatches, nLevels);
}

void PWCSolver::constructSe() {
	gf = greenOutside; // sets the global pointer to pass GF to C code
	apriori2_ = compression(&S_e_, waveletList, elementTree, nPatches, nLevels);
	WEM(&S_e_, waveletList, elementTree, T_, nPatches, nLevels, SingleLayer, DoubleLayer, -2*M_PI);
	aposteriori2_ = postproc(&S_e_, waveletList, elementTree, nPatches, nLevels);
}

void PWCSolver::solveFirstKind(const Eigen::VectorXd & potential, Eigen::VectorXd & charge) {
	double *rhs;
	double *u = (double*) calloc(nFunctions, sizeof(double));
    double * pot = const_cast<double *>(potential.data());
    double * chg = charge.data();
	double epsilon = greenOutside->getDielectricConstant();
	WEMRHS2M(&rhs, waveletList, elementTree, T_, nPatches, nLevels, pot,
			 quadratureLevel_);
	int iter = WEMPGMRES2(&S_i_, rhs, u, threshold, nPatches, nLevels);
	tdwtKon(u, nLevels, nFunctions);
	dwtKon(u, nLevels, nFunctions);
	for (int i = 0; i < nFunctions; i++) {
		rhs[i] += 4 * M_PI * u[i] / (epsilon - 1);
	}
	memset(u, 0, nFunctions * sizeof(double));
	iter = WEMPCG(&S_i_, rhs, u, threshold, nPatches, nLevels);
    tdwtKon(u, nLevels, nFunctions);
	energy_ext(u, pot, elementList, T_, nPatches, nLevels);
	charge_ext(u, chg, elementList, T_, nPatches, nLevels);
	free(rhs);
	free(u);
}

void PWCSolver::solveSecondKind(const Eigen::VectorXd & potential, Eigen::VectorXd & charge) 
{
	throw std::runtime_error("Second Kind not yet implemented"); //careful to the double layer sign when implementing it....
}

void PWCSolver::solveFull(const Eigen::VectorXd & potential, Eigen::VectorXd & charge) {
	double *rhs;
	double *u = (double*) calloc(nFunctions, sizeof(double));
	double *v = (double*) calloc(nFunctions, sizeof(double));
	//next line is just a quick fix to avoid problems with const but i do not like it...
    double * pot = const_cast<double *>(potential.data());
	WEMRHS2M(&rhs, waveletList, elementTree, T_, nPatches, nLevels, pot,
			 quadratureLevel_);
	int iters = WEMPCG(&S_i_, rhs, u, threshold, nPatches, nLevels);
	memset(rhs, 0, nFunctions*sizeof(double));
	for(unsigned int i = 0; i < nFunctions; i++) {
		for(unsigned int j = 0; j < S_e_.row_number[i]; j++)  {
			rhs[i] += S_e_.value1[i][j] * u[S_e_.index[i][j]];
		}
	}
	iters = WEMPGMRES3(&S_i_, &S_e_, rhs, v, threshold, 
					   nPatches, nLevels);
	for(unsigned int i = 0; i < nFunctions; i++) {
		u[i] -= 4*M_PI*v[i];
	}
	tdwtKon(u, nLevels, nFunctions);
	energy_ext(u, pot, elementList, T_, nPatches, nLevels);
	charge_ext(u, charge.data(), elementList, T_, nPatches, nLevels);
	free(rhs);
	free(u);
	free(v);
}

std::ostream & PWCSolver::printSolver(std::ostream & os) 
{
	os << "Solver Type: Wavelet, piecewise constant functions";
	return os;
}
