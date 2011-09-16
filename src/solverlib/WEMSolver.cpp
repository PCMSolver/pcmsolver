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
}

WEMSolver::WEMSolver(GreensFunction * gfi, GreensFunction * gfo) :
	PCMSolver(gfi, gfo) {
}

WEMSolver::WEMSolver(Section solver) : PCMSolver(solver) {
}

WEMSolver::~WEMSolver(){
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


