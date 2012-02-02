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

PCMSolver::PCMSolver(GreensFunctionInterface &gfi, GreensFunctionInterface &gfo){
	allocated = false;
	greenInside = &gfi; 
	greenOutside = &gfo;
}

PCMSolver::PCMSolver(GreensFunctionInterface *gfi, GreensFunctionInterface *gfo){
	allocated = false;
	greenInside = gfi; 
	greenOutside = gfo;
}

PCMSolver::PCMSolver(Section solver) {
	vector<Solvent> SolventData = PCMSolver::initSolvent();
	int solventIndex = solver.getInt("SolIndex");
	solvent = solver.getStr("Solvent");
	bool greenInsideSet  = solver.getSect("Green<inside>").isDefined();
	bool greenOutsideSet = solver.getSect("Green<outside>").isDefined();
	if (greenInsideSet) {
		greenInside = greenInside->allocateGreensFunctionInterface(solver.getSect("Green<inside>"));
	} else {
		greenInside = greenInside->allocateGreensFunctionInterface();
	}
	if (greenOutsideSet) {
		greenOutside = greenOutside->allocateGreensFunctionInterface(solver.getSect("Green<outside>"));
	} else {
		greenOutside = greenOutside->allocateGreensFunctionInterface(SolventData[solventIndex].getSolventEpsStatic());
	} 
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

vector<Solvent> PCMSolver::initSolvent() {
	/*
	  vector<Solvent> SolventData[] should contain all the solvent-related
	  data needed to set up the Green's functions and the non-electrostatic
	  terms calculations.
	  
	  These data are taken from the DALTON2011 internal implementation of
	  the Polarizable Continuum Model.
	*/

	vector<Solvent> SolventData(18);

	// ------------------------------------------------------------

	SolventData[0] = Solvent("Water", 78.39, 1.776, 1.385);
	SolventData[1] = Solvent("Methanol", 32.63, 1.758, 1.855);
	SolventData[2] = Solvent("Ethanol", 24.55, 1.847, 2.18);
	SolventData[3] = Solvent("Chloroform", 4.90, 2.085, 2.48);
	SolventData[4] = Solvent("Methylenechloride", 8.93, 2.020, 2.27);
	SolventData[5] = Solvent("1,2-Dichloroethane", 10.36, 2.085, 2.505);
	SolventData[6] = Solvent("Carbon tetrachloride", 2.228, 2.129, 2.685);
	SolventData[7] = Solvent("Benzene", 2.247, 2.244, 2.630);
	SolventData[8] = Solvent("Toluene", 2.379, 2.232, 2.82);
	SolventData[9] = Solvent("Chlorobenzene", 5.621, 2.320, 2.805);
	SolventData[10] = Solvent("Nitromethane", 38.20, 1.904, 2.155);
	SolventData[11] = Solvent("N-heptane", 1.92, 1.918, 3.125);
	SolventData[12] = Solvent("Cyclohexane", 2.023, 2.028, 2.815);
	SolventData[13] = Solvent("Aniline", 6.89, 2.506, 2.80);
	SolventData[14] = Solvent("Acetone", 20.7, 1.841, 2.38);
	SolventData[15] = Solvent("Tetrahydrofurane", 7.58, 1.971, 2.9);
	SolventData[16] = Solvent("Dimethylsulfoxide", 46.7, 2.179, 2.455);
	SolventData[17] = Solvent("Acetonitrile", 36.64, 1.806, 2.155);	
  
  // ------------------------------------------------------------

  return SolventData;
}

void PCMSolver::setSolvent(string & solv) {
	solvent = solv;
}

ostream & PCMSolver::printObject(ostream & os) {
	string type;
	if ( solverType == Traditional ) {
		type = "Traditional";
	} else {
		type = "Wavelet";
	}
	os << "~~~~~~~~~~ PCMSolver ~~~~~~~~~~\n" << endl;
	os << "========== Solver section" << endl;
	os << "Solver Type: " << type << endl;
	os << "Solvent: " << solvent << endl;
	return os;
}

