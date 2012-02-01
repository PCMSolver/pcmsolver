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
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "GreensFunctionSum.h"
#include "MetalSphere.h"
#include "Cavity.h"
#include "WaveletCavity.h"
#include "PCMSolver.h"
#include "WEMSolver.h"
#include "PWLSolver.h"

template <class T>
PWLSolver<T>::initPointers(GreensFunction<T> & gfi, GreensFunction<T> & gfo) : 
{
	elementTree = NULL;
	waveletList = NULL;
}

template <class T>
PWLSolver<T>::PWLSolver(GreensFunction<T> & gfi, GreensFunction<T> & gfo) : 
	WEMSolver<T>(gfi, gfo) {
	initPointers();
}

template <class T>
PWLSolver<T>::PWLSolver(GreensFunction<T> * gfi, GreensFunction<T> * gfo) :
	WEMSolver<T>(gfi, gfo) {
	initPointers();
}

template <class T>
PWLSolver<T>::PWLSolver(Section solver) : WEMSolver<T>(solver) {
	initPointers();
}

template <class T>
PWLSolver<T>::~PWLSolver(){
	if(elementTree != NULL) free_elementlist(&elementTree,nPatches,nLevels);
	if(waveletList != NULL) free_waveletlist(&waveletList,nPatches,nLevels);
}

template <class T>
void PWLSolver<T>::initInterpolation(WaveletCavity cavity) {
	init_interpolate(&T_, U, nPatches, nLevels);
	nNodes = gennet(&nodeList, &elementList, U, nPatches, nLevels);
}

template <class T>
void WEMSolver<T>::constructSystemMatrix(){
	generate_elementlist(&elementTree,nodeList,elementList,nPatches,nLevels);
	generate_waveletlist(&waveletList,elementTree,nPatches,nLevels);
	set_quadrature_level(waveletList,elementTree,nPatches,nLevels);
	simplify_waveletlist(waveletList,elementTree,nPatches,nLevels);
	complete_elementlist(waveletList,elementTree,nPatches,nLevels);
	
	this->fixPointersInside();
	apriori1_ = compression(&S_i_,waveletList,elementTree,nPatches,nLevels);
	WEM(&S_i_,waveletList,elementTree,T_,nPatches,nLevels,SLInt,DLUni,2*M_PI);
	//WEM(&S_i_,waveletList,elementTree,T_,nPatches,nLevels,SingleLayer,DoubleLayer,2*M_PI);
	aposteriori1_ = postproc(&S_i_,waveletList,elementTree,nPatches,nLevels);
	
	this->fixPointersOutside();
	apriori2_ = compression(&S_e_,waveletList,elementTree,nPatches,nLevels);
	WEM(&S_i_,waveletList,elementTree,T_,nPatches,nLevels,SLExt,DLUni,2*M_PI);
	//WEM(&S_e_,waveletList,elementTree,T_,nPatches,nLevels,SingleLayer,DoubleLayer,-2*M_PI);
	aposteriori2_ = postproc(&S_e_,waveletList,elementTree,nPatches,nLevels);
	systemMatricesInitialized_ = true;
}


template <class T>
void PWLSolver<T>::compCharge(const VectorXd & potential, VectorXd & charge) {
	double *rhs;
	double *u = (double*) calloc(nFunctions, sizeof(double));
	double *v = (double*) calloc(nFunctions, sizeof(double));
	//next line is just a quick fix but i do not like it...
    VectorXd pot = potential;
	WEMRHS2M(&rhs, waveletList, elementTree, T_, nPatches, nLevels, 
			 pot.data(), quadratureLevel_);
	int iters = WEMPCG(&S_i_, rhs, u, threshold, nPatches, nLevels);
	memset(rhs, 0, nFunctions*sizeof(double));
	for(unsigned int i = 0; i < nFunctions; i++) {
		for(unsigned int j = 0; j < S_e_.row_number[i]; j++)  {
			rhs[i] += S_e_.value1[i][j] * u[S_e_.index[i][j]];
		}
	}
	iters = WEMPGMRES3(&S_i_, &S_e_, rhs, v, threshold, nPatches, nLevels);
	for(unsigned int i=0; i<nFunctions; i++) {
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

template class PWLSolver <double>;
template class PWLSolver <taylor<double, 1, 1> >;
template class PWLSolver <taylor<double, 3, 1> >;
template class PWLSolver <taylor<double, 3 ,2> >;

