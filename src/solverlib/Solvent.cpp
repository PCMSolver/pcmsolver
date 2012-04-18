#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Solvent.h"

/*

  Methods for Solvent class
  written by Roberto Di Remigio, 2011

*/

Solvent::Solvent( const string & name, double epsStatic, 
				  double epsOptical, double radius ) {
	this->name = name;
	this->epsStatic = epsStatic;
	this->epsOptical = epsOptical;
	this->probeRadius = radius;
}

vector<Solvent> Solvent::initSolventVector() {
	/*
	  vector<Solvent> avaliableSolvents should contain all the solvent-related
	  data needed to set up the Green's functions and the non-electrostatic
	  terms calculations.
	  
	  These data are taken from the DALTON2011 internal implementation of
	  the Polarizable Continuum Model.
	*/

	vector<Solvent> solventData;
	solventData.push_back(Solvent("Water", 78.39, 1.776, 1.385));
	solventData.push_back(Solvent("Methanol", 32.63, 1.758, 1.855));
	solventData.push_back(Solvent("Ethanol", 24.55, 1.847, 2.18));
	solventData.push_back(Solvent("Chloroform", 4.90, 2.085, 2.48));
	solventData.push_back(Solvent("Methylenechloride", 8.93, 2.020, 2.27));
	solventData.push_back(Solvent("1,2-Dichloroethane", 10.36, 2.085, 2.505));
	solventData.push_back(Solvent("Carbon tetrachloride", 2.228, 2.129, 2.685));
	solventData.push_back(Solvent("Benzene", 2.247, 2.244, 2.630));
	solventData.push_back(Solvent("Toluene", 2.379, 2.232, 2.82));
	solventData.push_back(Solvent("Chlorobenzene", 5.621, 2.320, 2.805));
	solventData.push_back(Solvent("Nitromethane", 38.20, 1.904, 2.155));
	solventData.push_back(Solvent("N-heptane", 1.92, 1.918, 3.125));
	solventData.push_back(Solvent("Cyclohexane", 2.023, 2.028, 2.815));
	solventData.push_back(Solvent("Aniline", 6.89, 2.506, 2.80));
	solventData.push_back(Solvent("Acetone", 20.7, 1.841, 2.38));
	solventData.push_back(Solvent("Tetrahydrofurane", 7.58, 1.971, 2.9));
	solventData.push_back(Solvent("Dimethylsulfoxide", 46.7, 2.179, 2.455));
	solventData.push_back(Solvent("Acetonitrile", 36.64, 1.806, 2.155));	
	solventData.push_back(Solvent("Explicit", 0.0, 0.0, 0.0));	
  
  // ------------------------------------------------------------

  return solventData;
}

SolventMap Solvent::initSolventMap() {
	SolventMap mapSolvents;
	vector<Solvent> availSolvents = initSolventVector();
	for (int i = 0; i < availSolvents.size(); i++) {
		mapSolvents[availSolvents[i].getName()] = &availSolvents[i];
	}
	return mapSolvents;
}

ostream & Solvent::printObject(ostream & os) {
	string type = "Traditional";
	os << "Solvent name:           " << name << endl;
	os << "Static diel. constant:  " << epsStatic << endl;
	os << "Optical diel. constant: " << epsStatic << endl;
	os << "Solvent radius:         " << probeRadius;
	return os;
}


