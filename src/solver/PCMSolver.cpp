#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

#include "PCMSolver.h"

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
	return os;
}

