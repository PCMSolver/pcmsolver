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
#include "taylor.hpp"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "GreensFunctionSum.h"
#include "MetalSphere.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "PCMSolver.h"
#include "Solvent.h"

PCMSolver::PCMSolver(GreensFunctionInterface &gfi, GreensFunctionInterface &gfo,
					 int equation, int solver){
	allocated = false;
	greenInside = &gfi;
	greenOutside = &gfo;
	setSolverType(solver);
	setEquationType(equation);
}

PCMSolver::PCMSolver(GreensFunctionInterface *gfi, GreensFunctionInterface *gfo,
					 int equation, int solver){
	allocated = false;
	greenInside = gfi; 
	greenOutside = gfo;
	setSolverType(solver);
	setEquationType(equation);
}

PCMSolver::PCMSolver(const Section & solver) {
	SolventMap solvents = Solvent::initSolventMap();
	string name = solver.getStr("Solvent");
	this->solvent = *solvents["Water"];
	std::cout << this->solvent << std::endl;
	bool greenInsideSet  = solver.getSect("Green<inside>").isDefined();
	bool greenOutsideSet = solver.getSect("Green<outside>").isDefined();
	//	string PR = "ProbeRadius";
	//Keyword<double> radiusKeyword = solver.getKey(PR);
	//bool radiusSet = radiusKeyword.isDefined();
	if (greenInsideSet) {
		greenInside = greenInside->allocateGreensFunctionInterface(solver.getSect("Green<inside>"));
	} else {
		greenInside = greenInside->allocateGreensFunctionInterface();
	}
	if (greenOutsideSet) {
		greenOutside = greenOutside->allocateGreensFunctionInterface(solver.getSect("Green<outside>"));
	} else {
		greenOutside = greenOutside->allocateGreensFunctionInterface(this->solvent.getEpsStatic());
	}
	/*
	if (radiusSet) {
		this->solvent.setRadius(solver.getDbl("ProbeRadius"));
	}
	*/
	setSolverType(solver.getStr("SolverType"));
	setEquationType(solver.getStr("EquationType"));
	allocated = true;
}

PCMSolver::~PCMSolver(){
	if(allocated) {
		delete greenInside; 
		delete greenOutside;
	}
}

void PCMSolver::setSolverType(const string & type) {
	if (type == "Traditional") {
		setSolverType(Traditional);
	} else if (type == "Wavelet") {
		setSolverType(Wavelet);
	} else if (type == "Linear") {
		setSolverType(Linear);
	} else {
		exit(-1);
	}
}

void PCMSolver::setSolverType(int type) {
	switch (type) {
	case Traditional :
		solverType = Traditional;
		break;
	case Wavelet :
		solverType = Wavelet;
		break;
	case Linear :
		solverType = Linear;
		break;
	default : 
		exit(-1);
	}
}

void PCMSolver::setEquationType(const string & type) {
	if (type == "FirstKind") {
		setEquationType(FirstKind);
	} else if (type == "SecondKind") {
		setEquationType(SecondKind);
	} else if (type == "Full") {
		setEquationType(Full);
	} else {
		exit(-1);
	}
}

void PCMSolver::setEquationType(int type) {
	std::cout << "Setting eq type NUM" << type << std::endl; 
	switch (type) {
	case FirstKind :
		integralEquation = FirstKind;
		break;
	case SecondKind :
		integralEquation = SecondKind;
		break;
	case Full :
		integralEquation = Full;
		break;
	default : 
		exit(-1);
	}
}

GreensFunctionInterface & PCMSolver::getGreenInside(){
	return *greenInside;
}

GreensFunctionInterface & PCMSolver::getGreenOutside(){
	return *greenOutside;
}


GreensFunctionInterface * PCMSolver::getGreenInsideP(){
	return greenInside;
}

GreensFunctionInterface * PCMSolver::getGreenOutsideP(){
	return greenOutside;
}

void PCMSolver::setSolvent(Solvent & solv) {
	solvent = solv;
}

ostream & PCMSolver::printObject(ostream & os) {
	string type;
	if (solverType == Traditional) {
		type = "Traditional";
	} else if (solverType == Wavelet) {
		type = "Wavelets with piecewise constants";
	} else if (solverType == Linear) {
		type = "Wavelet with piecewise linears";
	} else {
		type = "Unknown";
	}
	os << "~~~~~~~~~~ PCMSolver ~~~~~~~~~~\n" << endl;
	os << "========== Solver section" << endl;
	os << "Solver Type: " << type << endl;
	os << "Solvent: " << solvent << endl;
	return os;
}

