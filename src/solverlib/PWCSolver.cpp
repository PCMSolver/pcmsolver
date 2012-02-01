/*! \file PWCSolver.cpp 
\brief PWC solver
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
#include "PWCSolver.h"

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

template <class T>
void PWCSolver<T>::initPointers()
{
	elementTree = NULL;
	waveletList = NULL;
}

template <class T>
PWCSolver<T>::PWCSolver(GreensFunction<T> & gfi, GreensFunction<T> & gfo) : 
	WEMSolver<T>(gfi, gfo) {
	initPointers();
}

template <class T>
PWCSolver<T>::PWCSolver(GreensFunction<T> * gfi, GreensFunction<T> * gfo) :
	WEMSolver<T>(gfi, gfo) {
	initPointers();
}

template <class T>
PWCSolver<T>::PWCSolver(Section solver) : WEMSolver<T>(solver) {
	initPointers();
}

template <class T>
PWCSolver<T>::~PWCSolver(){
	if(elementTree != NULL) free_elementlist(&elementTree,
											 this->nPatches,this->nLevels);
	if(waveletList != NULL) free_waveletlist(&waveletList,
											 this->nPatches,this->nLevels);
}

template <class T>
void PWCSolver<T>::initInterpolation() {
	init_interpolate(&this->T_, this->pointList, this->nPatches, this->nLevels);
	this->nNodes = gennet(&this->nodeList, &this->elementList, this->pointList, 
					this->nPatches, this->nLevels);
}

template <class T>
void PWCSolver<T>::constructSystemMatrix(){
	generate_elementlist(&elementTree,this->nodeList,this->elementList,
						 this->nPatches, this->nLevels);
	generate_waveletlist(&waveletList,elementTree,this->nPatches,
						 this->nLevels);
	set_quadrature_level(waveletList,elementTree,this->nPatches,
						 this->nLevels);
	simplify_waveletlist(waveletList,elementTree,this->nPatches,
						 this->nLevels);
	complete_elementlist(waveletList,elementTree,this->nPatches,
						 this->nLevels);
	
	this->fixPointersInside();
	this->apriori1_ = compression(&this->S_i_,waveletList,elementTree,
							this->nPatches,this->nLevels);
	WEM(&this->S_i_,waveletList,elementTree,this->T_,this->nPatches,
		this->nLevels,SLInt,DLUni,2*M_PI);
	this->aposteriori1_ = postproc(&this->S_i_,waveletList,elementTree,
							 this->nPatches,this->nLevels);
	
	this->fixPointersOutside();
	this->apriori2_ = compression(&this->S_e_,waveletList,elementTree,
								  this->nPatches,this->nLevels);
	WEM(&this->S_i_,waveletList,elementTree,this->T_,this->nPatches,
		this->nLevels,SLExt,DLUni,2*M_PI);
	this->aposteriori2_ = postproc(&this->S_e_,waveletList,elementTree,
							 this->nPatches,this->nLevels);
	this->systemMatricesInitialized_ = true;
}


template <class T>
void PWCSolver<T>::compCharge(const VectorXd & potential, VectorXd & charge) {
	double *rhs;
	double *u = (double*) calloc(this->nFunctions, sizeof(double));
	double *v = (double*) calloc(this->nFunctions, sizeof(double));
	//next line is just a quick fix but i do not like it...
    VectorXd pot = potential;
	WEMRHS2M(&rhs, waveletList, elementTree, this->T_, this->nPatches, this->nLevels, 
			 pot.data(), this->quadratureLevel_);
	int iters = WEMPCG(&this->S_i_, rhs, u, this->threshold, this->nPatches, this->nLevels);
	memset(rhs, 0, this->nFunctions*sizeof(double));
	for(unsigned int i = 0; i < this->nFunctions; i++) {
		for(unsigned int j = 0; j < this->S_e_.row_number[i]; j++)  {
			rhs[i] += this->S_e_.value1[i][j] * u[this->S_e_.index[i][j]];
		}
	}
	iters = WEMPGMRES3(&this->S_i_, &this->S_e_, rhs, v, this->threshold, 
					   this->nPatches, this->nLevels);
	for(unsigned int i = 0; i < this->nFunctions; i++) {
		u[i] -= 4*M_PI*v[i];
	}
	tdwtKon(u, this->nLevels, this->nFunctions);
  // Interpolate charges
	cubature *Q;
	init_Gauss_Square(&Q, this->quadratureLevel_+1);
	vector2 s;
	vector2 t;
	int index = 0;
	int zi = 0;
	double h = 1.0/(1<<this->nLevels);
	for (unsigned int i1=0; i1<this->nPatches; i1++) {
	    s.y = 0;
		for (int i2=0; i2<(1<<this->nLevels); i2++) {
			s.x = 0;
			for (int i3=0; i3<(1<<this->nLevels); i3++) {
				for (unsigned int k=0; k<Q[this->quadratureLevel_].nop; k++) {
					t = vector2_add(s, vector2_Smul(h, Q[this->quadratureLevel_].xi[k]));
					charge(index) = Q[this->quadratureLevel_].w[k]*u[zi] * h;
					index ++;
				}
				s.x += h;
				zi++;
			}
			s.y += h;
		}
	}
	free_Gauss_Square(&Q, this->quadratureLevel_ + 1);
	free(rhs);
	free(u);
	free(v);
	charge /= -ToAngstrom; //WARNING  WARNING  WARNING
}

template class PWCSolver <double>;
template class PWCSolver <taylor<double, 1, 1> >;
template class PWCSolver <taylor<double, 3, 1> >;
template class PWCSolver <taylor<double, 3 ,2> >;

