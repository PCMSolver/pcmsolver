#include "WEMSolver.hpp"

#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif
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
}

#include "Cavity.hpp"
#include "GreensFunction.hpp"
#include "WaveletCavity.hpp"

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

/*
WEMSolver::WEMSolver(const Section & solver) : PCMSolver(solver) {
	initWEMMembers();
}*/

WEMSolver::~WEMSolver()
{
	if(nodeList != NULL)    free(nodeList);
	if(elementList != NULL) free_patchlist(&elementList,nFunctions);
	if(pointList != NULL)   free_points(&pointList, nPatches, nLevels);
	if(systemMatricesInitialized_){
		free_sparse2(&S_i_);
		if(integralEquation == Full) free_sparse2(&S_e_);
	}
}

void WEMSolver::uploadCavity(WaveletCavity & cavity) {
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
				Eigen::Vector3d p = cavity.getNodePoint(kk);
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
		this->constructWavelets();
		this->constructSystemMatrix();
	} else {
		throw std::runtime_error("Wavelet type cavity needed for wavelet solver.");
	}
}

void WEMSolver::constructSystemMatrix(){
	constructSi();
	if(integralEquation == Full) {
		constructSe();
	}
}

void WEMSolver::compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge) {

	switch (integralEquation) {
	case FirstKind:
		solveFirstKind(potential, charge);
		break;
	case SecondKind:
		solveSecondKind(potential, charge);
		break;
	case Full:
		solveFull(potential, charge);
		break;
	default:
		throw std::runtime_error("Invalid case");
	}
	charge *= -1.0;
	//	charge /= -ToAngstrom; //WARNING  WARNING  WARNING
}
