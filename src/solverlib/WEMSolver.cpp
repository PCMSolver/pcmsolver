/*! \file WEMSolver.cpp 
\brief WEM solver
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
}

#include "Getkw.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "GreensFunctionSum.h"
#include "MetalSphere.h"
#include "Cavity.h"
#include "WaveletCavity.h"
#include "PCMSolver.h"
#include "WEMSolver.h"

static GreensFunction *gf = NULL;

static double SingleLayer(vector3 x, vector3 y){
  
  Vector3d vx(x.x, x.y, x.z);
  Vector3d vy(y.x, y.y, y.z);
  
  return gf->evalf(vx, vy);
}

static double DoubleLayer(vector3 x, vector3 y, vector3 n_y){
  
  Vector3d vx(x.x, x.y, x.z);
  Vector3d vy(y.x, y.y, y.z);
  Vector3d vn_y(n_y.x, n_y.y, n_y.z);
  
  return gf->evald(vn_y, vx, vy);
}
WEMSolver::WEMSolver(GreensFunction & gfi, GreensFunction & gfo) : 
	PCMSolver(gfi, gfo) {
	nodeList = NULL;
	elementList = NULL;
	T_ = NULL;
	elementTree = NULL;
	waveletList = NULL;
	systemMatricesInitialized_ = false;
	threshold = 1e-10;
	quadratureLevel_ = 1;
	nQuadPoints = 0;
}

WEMSolver::WEMSolver(GreensFunction * gfi, GreensFunction * gfo) :
	PCMSolver(gfi, gfo) {
	nodeList = NULL;
	elementList = NULL;
	T_ = NULL;
	elementTree = NULL;
	waveletList = NULL;
	systemMatricesInitialized_ = false;
	threshold = 1e-10;
	quadratureLevel_ = 1;
	nQuadPoints = 0;
}

WEMSolver::WEMSolver(Section solver) : PCMSolver(solver) {
	nodeList = NULL;
	elementList = NULL;
	T_ = NULL;
	elementTree = NULL;
	waveletList = NULL;
	systemMatricesInitialized_ = false;
	threshold = 1e-10;
	quadratureLevel_ = 1;
	nQuadPoints = 0;
}

WEMSolver::~WEMSolver(){
	if(nodeList != NULL) free(nodeList);
	if(elementList != NULL) free_patchlist(&elementList,nFunctions);
	if(T_ != NULL) free_interpolate(&T_,nPatches,nLevels);
	if(elementTree != NULL) free_elementlist(&elementTree,nPatches,nLevels);
	if(waveletList != NULL) free_waveletlist(&waveletList,nPatches,nLevels);
	
	if(systemMatricesInitialized_){
		free_sparse2(&S_i_);
		free_sparse2(&S_e_);
	}
}

void WEMSolver::uploadCavity(WaveletCavity cavity) {
	
	vector3 ***U = NULL;
	nPatches = cavity.getNPatches();
	nLevels = cavity.getNLevels();

	cout << "nPatches " << nPatches << endl;
	cout << "nLevels " << nLevels << endl;
	int n = 1<<nLevels;
	nFunctions = nPatches * n * n;
	alloc_points(&U, nPatches, nLevels);
	int kk = 0;
	cout << "Input init interpolate " << nPatches << " " << nLevels <<endl;
	for (int i = 0; i < nPatches; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= n; k++) {
				Vector3d p = cavity.getNodePoint(kk);
				U[i][k][j] = vector3_make(p(0), p(1), p(2));
				kk++;
			}
		}
	}	
	for (int i = 0; i < nPatches; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k < n; k++) {
				std::cout << " " << i << " " << " " << j << " " << k << " ";
				std::cout << U[i][j][k].x << " " << U[i][j][k].y << " " << U[i][j][k].z << std::endl;
			}
		}
	}
	init_interpolate(&T_, U, nPatches, nLevels);
	nNodes = gennet(&nodeList, &elementList, U, nPatches, nLevels);
	free_points(&U, nPatches, nLevels);
}

void WEMSolver::constructSystemMatrix(){

  generate_elementlist(&elementTree,nodeList,elementList,nPatches,nLevels);
  generate_waveletlist(&waveletList,elementTree,nPatches,nLevels);
  set_quadrature_level(waveletList,elementTree,nPatches,nLevels);
  simplify_waveletlist(waveletList,elementTree,nPatches,nLevels);
  complete_elementlist(waveletList,elementTree,nPatches,nLevels);
  
  // Default GFs ok for the internal part
  gf = getGreenInsideP();
  apriori1_ = compression(&S_i_,waveletList,elementTree,nPatches,nLevels);
  WEM(&S_i_,waveletList,elementTree,T_,nPatches,nLevels,SingleLayer,DoubleLayer,2*M_PI);
  aposteriori1_ = postproc(&S_i_,waveletList,elementTree,nPatches,nLevels);

  // Use interface for the external part
  gf = getGreenOutsideP();
  apriori2_ = compression(&S_e_,waveletList,elementTree,nPatches,nLevels);
  WEM(&S_e_,waveletList,elementTree,T_,nPatches,nLevels,SingleLayer,DoubleLayer,-2*M_PI);
  aposteriori2_ = postproc(&S_e_,waveletList,elementTree,nPatches,nLevels);
  
  systemMatricesInitialized_ = true;
}


VectorXd WEMSolver::compCharge(const VectorXd &potential) {
	VectorXd charge;
	cout << "WEM solver not yet implemented" << endl;
	exit(1);
	return charge;
}

void WEMSolver::compCharge(const VectorXd & potential, VectorXd & charge) {
	double *rhs;
	double *u = (double*) calloc(nFunctions, sizeof(double));
	double *v = (double*) calloc(nFunctions, sizeof(double));
	
	//next line is just a quick fix but i do not like it...
    VectorXd pot = potential;

	cout << "comp charges 1" << endl;

	cout << " " << nPatches << " " << nLevels << " " << quadratureLevel_ << endl;
	cout << pot.transpose() << endl;

	WEMRHS2M(&rhs, waveletList, elementTree, T_, nPatches, nLevels, 
			 pot.data(), quadratureLevel_);

	cout << "comp charges 2" << endl;

	int iters = WEMPCG(&S_i_, rhs, u, threshold, nPatches, nLevels);
	cout << "comp charges 3" << endl;

	memset(rhs, 0, nFunctions*sizeof(double));
	cout << "comp charges 4" << endl;
	for(unsigned int i = 0; i < nFunctions; i++) {
		cout << i << " " << S_e_.row_number[i] << " " << nFunctions << endl; 
		for(unsigned int j = 0; j < S_e_.row_number[i]; j++)  {
			rhs[i] += S_e_.value1[i][j] * u[S_e_.index[i][j]];
			cout << "comp charges ij " << i << " " << j << " " << rhs[i] << endl;
		}
	}


	cout << "comp charges 5" << endl;
	iters = WEMPGMRES3(&S_i_, &S_e_, rhs, v, threshold, nPatches, nLevels);
	cout << "comp charges 6" << endl;
	for(unsigned int i=0; i<nFunctions; i++) u[i] -= 4*M_PI*v[i];
	cout << "comp charges 7" << endl;
	tdwtKon(u, nLevels, nFunctions);
	cout << "comp charges 8" << endl;

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
}
