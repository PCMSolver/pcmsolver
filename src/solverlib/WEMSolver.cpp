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

static double (*SingleLayer) (vector3 x, vector3 y);
static double (*DoubleLayer) (vector3 x, vector3 y, vector3 n_y);

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


typedef taylor<double, 3, 1> T1;
typedef taylor<double, 3, 2> T2;
typedef taylor<double, 1, 1> T3;

static WEMSolver<double> * globalSolverD;
static WEMSolver<T1> * globalSolverT1;
static WEMSolver<T2> * globalSolverT2;
static WEMSolver<T3> * globalSolverT3;


static double SingleLayerD (vector3 x, vector3 y) {
	return globalSolverD->SL(x, y);
}
static double DoubleLayerD (vector3 x, vector3 y, vector3 n_y) {
	double val = globalSolverD->DL(x, y, n_y);
	return val;
}
static double SingleLayerT1 (vector3 x, vector3 y) {
	return globalSolverT1->SL(x, y);
}
static double DoubleLayerT1 (vector3 x, vector3 y, vector3 n_y) {
	double val = globalSolverT1->DL(x, y, n_y);
	return val;
}
static double SingleLayerT2 (vector3 x, vector3 y) {
	return globalSolverT2->SL(x, y);
}
static double DoubleLayerT2 (vector3 x, vector3 y, vector3 n_y) {
	double val = globalSolverT2->DL(x, y, n_y);
	return val;
}
static double SingleLayerT3 (vector3 x, vector3 y) {
	return globalSolverT3->SL(x, y);
}
static double DoubleLayerT3 (vector3 x, vector3 y, vector3 n_y) {
	double val = globalSolverT3->DL(x, y, n_y);
	return val;
}

template <class T>
void WEMSolver<T>::fixPointersInside() {
	this->gf = this->greenInside;
	if (typeid(this) == typeid(WEMSolver<double> *)) {
		globalSolverD = dynamic_cast<WEMSolver<double> * > (this);
		SingleLayer = &SingleLayerD;
		DoubleLayer = &DoubleLayerD;
	} else if (typeid(this) == typeid(WEMSolver<T1> *)) {
		globalSolverT1 = dynamic_cast<WEMSolver<T1> * > (this);
		SingleLayer = &SingleLayerT1;
		DoubleLayer = &DoubleLayerT1;
    } else if (typeid(this) == typeid(WEMSolver<T2> *)) {
		globalSolverT2 = dynamic_cast<WEMSolver<T2> * > (this);
		SingleLayer = &SingleLayerT2;
		DoubleLayer = &DoubleLayerT2;
    } else if (typeid(this) == typeid(WEMSolver<T3> *)) {
		globalSolverT3 = dynamic_cast<WEMSolver<T3> * > (this);
		SingleLayer = &SingleLayerT3;
		DoubleLayer = &DoubleLayerT3;
	} else {
		std::cout << "this is not good! "<< std::endl;
		exit(-1);
	}
}

template <class T>
void WEMSolver<T>::fixPointersOutside() {
	this->gf = this->greenOutside;
	if (typeid(this) == typeid(WEMSolver<double> *)) {
		globalSolverD = dynamic_cast<WEMSolver<double> * > (this);
		SingleLayer = &SingleLayerD;
		DoubleLayer = &DoubleLayerD;
	} else if (typeid(this) == typeid(WEMSolver<T1> *)) {
		globalSolverT1 = dynamic_cast<WEMSolver<T1> * > (this);
		SingleLayer = &SingleLayerT1;
		DoubleLayer = &DoubleLayerT1;
	} else if (typeid(this) == typeid(WEMSolver<T2> *)) {
		globalSolverT2 = dynamic_cast<WEMSolver<T2> * > (this);
		SingleLayer = &SingleLayerT2;
		DoubleLayer = &DoubleLayerT2;
	} else if (typeid(this) == typeid(WEMSolver<T3> *)) {
		globalSolverT3 = dynamic_cast<WEMSolver<T3> * > (this);
		SingleLayer = &SingleLayerT3;
		DoubleLayer = &DoubleLayerT3;
	} else {
		std::cout << "this is not good! "<< std::endl;
		exit(-1);
	}
}

template <class T>
double WEMSolver<T>::SL(vector3 x, vector3 y){
  Vector3d vx(x.x, x.y, x.z);
  Vector3d vy(y.x, y.y, y.z);
  double value = this->gf->evalf(vx, vy);
  return value;
}

template <class T>
double WEMSolver<T>::DL(vector3 x, vector3 y, vector3 n_y){  
  Vector3d vx(x.x, x.y, x.z);
  Vector3d vy(y.x, y.y, y.z);
  Vector3d vn_y(n_y.x, n_y.y, n_y.z);
  double value = this->gf->evald(vn_y, vx, vy);
  return value;
}

template <class T>
void WEMSolver<T>::initWEMMembers()
{
	pointList = NULL;
	nodeList = NULL;
	elementList = NULL;
	T_ = NULL;
	systemMatricesInitialized_ = false;
	threshold = 1e-10;
	quadratureLevel_ = 1;
	nQuadPoints = 0;
}

template <class T>
WEMSolver<T>::WEMSolver(GreensFunction<T> & gfi, GreensFunction<T> & gfo) : 
	PCMSolver<T>(gfi, gfo) {
	initWEMMembers();
}

template <class T>
WEMSolver<T>::WEMSolver(GreensFunction<T> * gfi, GreensFunction<T> * gfo) :
	PCMSolver<T>(gfi, gfo) {
	initWEMMembers();
}

template <class T>
WEMSolver<T>::WEMSolver(Section solver) : PCMSolver<T>(solver) {
	initWEMMembers();
}

template <class T>
WEMSolver<T>::~WEMSolver(){
	if(nodeList != NULL)    free(nodeList);
	if(elementList != NULL) free_patchlist(&elementList,nFunctions);
	if(T_ != NULL)          free_interpolate(&T_,nPatches,nLevels);
	if(pointList != NULL)   free_points(&pointList, nPatches, nLevels);
	if(systemMatricesInitialized_){
		free_sparse2(&S_i_);
		free_sparse2(&S_e_);
	}
}

template <class T>
void WEMSolver<T>::uploadCavity(WaveletCavity cavity) {
	nPatches = cavity.getNPatches();
	nLevels = cavity.getNLevels();
	int n = (1<<nLevels);
	nFunctions = nPatches * n * n;
	alloc_points(&pointList, nPatches, nLevels);
	int kk = 0;
	// Ask Helmut about index switch
	for (int i = 0; i < nPatches; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k <= n; k++) {
				Vector3d p = cavity.getNodePoint(kk);
				pointList[i][k][j] = vector3_make(p(0), p(1), p(2));
				kk++;
			}
		}
	}	
}

template <class T>
void WEMSolver<T>::buildSystemMatrix(Cavity & cavity) {
    if (WaveletCavity *waveletCavity = dynamic_cast<WaveletCavity*> (&cavity)) {
		this->uploadCavity(*waveletCavity);
		this->initInterpolation();
		this->constructSystemMatrix();
	} else {
		std::cout << "Wavelet-type cavity needed for wavelet solver." 
				  << std::endl;
		exit(-1);
	}
}

template class WEMSolver <double>;
template class WEMSolver <taylor<double, 1, 1> >;
template class WEMSolver <taylor<double, 3, 1> >;
template class WEMSolver <taylor<double, 3 ,2> >;

