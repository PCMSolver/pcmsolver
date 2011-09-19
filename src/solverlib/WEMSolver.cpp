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

#include "Getkw.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "GreensFunctionSum.h"
#include "MetalSphere.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "PCMSolver.h"
#include "WEMSolver.h"

WEMSolver::WEMSolver(GreensFunction & gfi, GreensFunction & gfo) : 
	PCMSolver(gfi, gfo) {
	P_ = NULL;
	F_ = NULL;
	T_ = NULL;
	E_ = NULL;
	W_ = NULL;
	systemMatricesInitialized_ = false;
	eps_ = 1e-10;
	quadratureLevel_ = 1;
	nPoints_ = 0;
}

WEMSolver::WEMSolver(GreensFunction * gfi, GreensFunction * gfo) :
	PCMSolver(gfi, gfo) {
	P_ = NULL;
	F_ = NULL;
	T_ = NULL;
	E_ = NULL;
	W_ = NULL;
	systemMatricesInitialized_ = false;
	eps_ = 1e-10;
	quadratureLevel_ = 1;
	nPoints_ = 0;
}

WEMSolver::WEMSolver(Section solver) : PCMSolver(solver) {
	P_ = NULL;
	F_ = NULL;
	T_ = NULL;
	E_ = NULL;
	W_ = NULL;
	systemMatricesInitialized_ = false;
	eps_ = 1e-10;
	quadratureLevel_ = 1;
	nPoints_ = 0;
}

WEMSolver::~WEMSolver(){
	if(P_ != NULL) free(P_);
	if(F_ != NULL) free_patchlist(&F_,nf_);
	if(T_ != NULL) free_interpolate(&T_,p_,M_);
	if(E_ != NULL) free_elementlist(&E_,p_,M_);
	if(W_ != NULL) free_waveletlist(&W_,p_,M_);
	
	if(systemMatricesInitialized_){
		free_sparse2(&S_i_);
		free_sparse2(&S_e_);
	}
}

VectorXd WEMSolver::compCharge(const VectorXd &potential) {
	VectorXd charge;
	cout << "WEM solver not yet implemented" << endl;
	exit(1);
	return charge;
}

void WEMSolver::compCharge(const VectorXd & potential, VectorXd & charge) {
	cout << "WEM solver not yet implemented" << endl;
	exit(1);
}


