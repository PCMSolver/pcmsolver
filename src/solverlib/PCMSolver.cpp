/*! \file IEFSolver.cpp 
\brief PCM solver
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

PCMSolver::PCMSolver(GreensFunction &gfi, GreensFunction &gfo){
	allocated = false;
	greenInside = &gfi; 
	greenOutside = &gfo;
}

PCMSolver::PCMSolver(GreensFunction *gfi, GreensFunction *gfo){
	allocated = false;
	greenInside = gfi; 
	greenOutside = gfo;
}

PCMSolver::PCMSolver(Section solver) {
	allocated = true;
	greenInside  = 
		greenInside->allocateGreensFunction(solver.getSect("Green<inside>"));
	greenOutside = 
		greenOutside->allocateGreensFunction(solver.getSect("Green<outside>"));
}

PCMSolver::~PCMSolver(){
	if(allocated) {
		delete greenInside; 
		delete greenOutside;
	}
}

void PCMSolver::setSolverType(const string & type) {
	if (type == "Traditional") {
		this->setSolverType(Traditional);
	} else if (type == "Wavelet") {
		this->setSolverType(Wavelet);
	} else {
		exit(-1);
	}
}

void PCMSolver::setSolverType(int type) {
	switch (type) {
	case Traditional :
		this->solverType = Traditional;
		break;
	case Wavelet :
		this->solverType = Wavelet;
		break;
	default : 
		exit(-1);
	}
}

GreensFunction & PCMSolver::getGreenInside(){
	return *greenInside;
}

GreensFunction & PCMSolver::getGreenOutside(){
	return *greenOutside;
}


GreensFunction * PCMSolver::getGreenInsideP(){
	return greenInside;
}

GreensFunction * PCMSolver::getGreenOutsideP(){
	return greenOutside;
}

