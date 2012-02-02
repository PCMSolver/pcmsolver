/*! \file PWLSolver.cpp 
\brief PWL solver
*/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

extern "C"{
#include "vector3.h"
#include "sparse2.h"
#include "intvector_pwl.h"
#include "basis_pwl.h"
#include "WEM_pwl.h"
#include "read_points.h"
#include "vector2.h"
#include "interpolate_pwl.h"
#include "topology_pwl.h"
#include "kern.h"
#include "compression_pwl.h"
#include "postproc_pwl.h"
#include "WEMRHS_pwl.h"
#include "WEMPCG_pwl.h"
#include "WEMPGMRES_pwl.h"
#include "dwt_pwl.h"
#include "cubature.h"
#include "gauss_square.h"
#include "constants.h"
}

#include "Constants.h"
#include "Getkw.h"
#include "taylor.hpp"
#include "GreensFunctionInterface.h"
#include "Cavity.h"
#include "WaveletCavity.h"
#include "PCMSolver.h"
#include "WEMSolver.h"
#include "PWLSolver.h"

static GreensFunctionInterface * gf;

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
	Vector3d grad, dir;
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

static double SingleLayer (vector3 x, vector3 y) {
	Vector3d vx(x.x, x.y, x.z);
	Vector3d vy(y.x, y.y, y.z);
	double value = gf->evalf(vx, vy);
	return value;
}

static double DoubleLayer (vector3 x, vector3 y, vector3 n_y) {
	Vector3d vx(x.x, x.y, x.z);
	Vector3d vy(y.x, y.y, y.z);
	Vector3d vn_y(n_y.x, n_y.y, n_y.z);
	double value = gf->evald(vn_y, vx, vy);
	return value;
}

void PWLSolver::initPointers()
{
	elementTree = NULL;
	waveletList = NULL;
}

PWLSolver::PWLSolver(GreensFunctionInterface & gfi, GreensFunctionInterface & gfo) : 
	WEMSolver(gfi, gfo) {
	initPointers();
}

PWLSolver::PWLSolver(GreensFunctionInterface * gfi, GreensFunctionInterface * gfo) :
	WEMSolver(gfi, gfo) {
	initPointers();
}

PWLSolver::PWLSolver(Section solver) : WEMSolver(solver) {
	initPointers();
}

PWLSolver::~PWLSolver(){
	if(elementTree != NULL) free_elementlist_pwl(&elementTree, nPatches,nLevels);
	if(waveletList != NULL) free_waveletlist_pwl(&waveletList, nPatches,nLevels);
}

void PWLSolver::initInterpolation() {
	init_interpolate(&T_, pointList, nPatches, nLevels);
	nNodes = gennet(&nodeList, &elementList, pointList, 
					nPatches, nLevels);
}

void PWLSolver::constructSystemMatrix(){
	generate_elementlist_pwl(&elementTree, nodeList, elementList, nPatches, nLevels);
	generate_waveletlist_pwl(&waveletList, elementTree, nPatches, nLevels);
	set_quadrature_level_pwl(waveletList,elementTree,nPatches, nLevels);
	simplify_waveletlist_pwl(waveletList,elementTree,nPatches, nLevels);
	complete_elementlist_pwl(waveletList,elementTree,nPatches, nLevels);
	
	gf = greenInside; // sets the global pointer to pass GF to C code
	apriori1_ = compression_pwl(&S_i_,waveletList,elementTree, nPatches,nLevels);
	WEM_pwl(&S_i_,waveletList,elementTree,T_,nPatches, nLevels,SLInt,DLUni,2*M_PI);
	aposteriori1_ = postproc_pwl(&S_i_,waveletList,elementTree, nPatches,nLevels);
	
	gf = greenOutside; // sets the global pointer to pass GF to C code
	apriori2_ = compression_pwl(&S_e_,waveletList,elementTree, nPatches,nLevels);
	WEM_pwl(&S_e_,waveletList,elementTree,T_,nPatches,nLevels,SLExt,DLUni,-2*M_PI);
	aposteriori2_ = postproc_pwl(&S_e_,waveletList,elementTree, nPatches,nLevels);
	systemMatricesInitialized_ = true;
}

void PWLSolver::compCharge(const VectorXd & potential, VectorXd & charge) {
	double *rhs;
	double *u = (double*) calloc(nFunctions, sizeof(double));
	double *v = (double*) calloc(nFunctions, sizeof(double));
	//next line is just a quick fix but i do not like it...
    VectorXd pot = potential;
	WEMRHS2M_pwl(&rhs, waveletList, elementTree, T_, nPatches, nLevels, 
			 pot.data(), quadratureLevel_);
	int iters = WEMPCG_pwl(&S_i_, rhs, u, threshold, nPatches, nLevels);
	memset(rhs, 0, nFunctions*sizeof(double));
	for(unsigned int i = 0; i < nFunctions; i++) {
		for(unsigned int j = 0; j < S_e_.row_number[i]; j++)  {
			rhs[i] += S_e_.value1[i][j] * u[S_e_.index[i][j]];
		}
	}
	iters = WEMPGMRES3_pwl(&S_i_, &S_e_, rhs, v, threshold, 
					   nPatches, nLevels);
	for(unsigned int i = 0; i < nFunctions; i++) {
		u[i] -= 4*M_PI*v[i];
	}
	tdwtKon(u, nLevels, nFunctions);
  // Interpolate charges
	cubature *Q;
	init_Gauss_Square(&Q, quadratureLevel_+1);
	vector2 s;
	vector2 t;
	int index = 0;
	int zi = 0;
	double h = 1.0/(1<<nLevels);
	for (unsigned int i1=0; i1<nPatches; i1++) {
	    s.y = 0;
		for (int i2=0; i2<(1<<nLevels); i2++) {
			s.x = 0;
			for (int i3=0; i3<(1<<nLevels); i3++) {
				for (unsigned int k=0; k<Q[quadratureLevel_].nop; k++) {
					t = vector2_add(s, vector2_Smul(h, Q[quadratureLevel_].xi[k]));
					charge(index) = Q[quadratureLevel_].w[k]*u[zi] * h;
					index ++;
				}
				s.x += h;
				zi++;
			}
			s.y += h;
		}
	}
	free_Gauss_Square(&Q, quadratureLevel_ + 1);
	free(rhs);
	free(u);
	free(v);
	charge /= -ToAngstrom; //WARNING  WARNING  WARNING
}

