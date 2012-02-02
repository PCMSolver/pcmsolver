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
#include "GreensFunctionInterface.h"
#include "Cavity.h"
#include "WaveletCavity.h"
#include "PCMSolver.h"
#include "WEMSolver.h"

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
	c.x = x.x-y.x;
	c.y = x.y-y.y;
	c.z = x.z-y.z;
	r = sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
	return (c.x*n_y.x+c.y*n_y.y+c.z*n_y.z)/(r*r*r);
}

void WEMSolver::initWEMMembers()
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

WEMSolver::WEMSolver(GreensFunctionInterface & gfi, GreensFunctionInterface & gfo) : 
	PCMSolver(gfi, gfo) {
	initWEMMembers();
}

WEMSolver::WEMSolver(GreensFunctionInterface * gfi, GreensFunctionInterface * gfo) :
	PCMSolver(gfi, gfo) {
	initWEMMembers();
}

WEMSolver::WEMSolver(Section solver) : PCMSolver(solver) {
	initWEMMembers();
}

WEMSolver::~WEMSolver(){
	if(nodeList != NULL)    free(nodeList);
	if(elementList != NULL) free_patchlist(&elementList,nFunctions);
	if(T_ != NULL)          free_interpolate(&T_,nPatches,nLevels);
	if(pointList != NULL)   free_points(&pointList, nPatches, nLevels);
	if(systemMatricesInitialized_){
		free_sparse2(&S_i_);
		free_sparse2(&S_e_);
	}
}

void WEMSolver::uploadCavity(WaveletCavity cavity) {
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

void WEMSolver::buildSystemMatrix(Cavity & cavity) {
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

