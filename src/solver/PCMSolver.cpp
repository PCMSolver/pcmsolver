#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

//#include "Getkw.h"
//#include "taylor.hpp"
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

/*
PCMSolver::PCMSolver(const Section & solver) {
	std::map<std::string, Solvent> solvents = Solvent::initSolventMap();
	string name = solver.getStr("Solvent");
	this->solvent = solvents[name];
        if (name == "Explicit") {
		greenInside = greenInside->allocateGreensFunctionInterface(solver.getSect("Green<inside>"));
	} else {
		greenInside = greenInside->allocateGreensFunctionInterface();
	}
        if (name == "Explicit") {
		greenOutside = greenOutside->allocateGreensFunctionInterface(solver.getSect("Green<outside>"));
	} else {
		greenOutside = greenOutside->allocateGreensFunctionInterface(this->solvent.getEpsStatic());
	}
	setSolverType(solver.getStr("SolverType"));
	setEquationType(solver.getStr("EquationType"));
	allocated = true;
}*/

PCMSolver::~PCMSolver(){
	if(allocated) {
		delete greenInside; 
		delete greenOutside;
	}
}

void PCMSolver::setSolverType(const std::string & type) {
	if (type == "IEFPCM") {
		setSolverType(IEFPCM);
	} else if (type == "CPCM") {
		setSolverType(CPCM);
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
	case IEFPCM :
		solverType = IEFPCM;
		break;
	case CPCM :
		solverType = CPCM;
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

void PCMSolver::setEquationType(const std::string & type) {
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
	switch (type) {
	case FirstKind :
		integralEquation = FirstKind;
		break;
	case SecondKind :
		std::cout << "Setting eq type NUM " << type << std::endl; 
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

std::ostream & PCMSolver::printObject(std::ostream & os) {
	std::string type;
	if (solverType == IEFPCM) {
		type = "IEFPCM";
	} else if (solverType == CPCM) {
		type = "CPCM";
	} else if (solverType == Wavelet) {
		type = "Wavelets with piecewise constants";
	} else if (solverType == Linear) {
		type = "Wavelet with piecewise linears";
	} else {
		type = "Unknown";
	}
	os << "~~~~~~~~~~ PCMSolver ~~~~~~~~~~\n" << std::endl;
	os << "========== Solver section" << std::endl;
	os << "Solver Type: " << type << std::endl;
	os << solvent; // Solvent has its own operator<< overloaded
	return os;
}

