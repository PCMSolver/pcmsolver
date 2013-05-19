#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "SurfaceFunction.h"
#include "Cavity.h"

/*

Methods for basic cavity class
written by Krzysztof Mozgawa, 2011

*/

void Cavity::writeOutput(string &filename){
    ofstream output;
    output.open(filename.c_str(), fstream::out);
    output << nTess << endl;
    for(int i=0; i < nTess; i++) {
		output << tessCenter(0,i) << " ";
		output << tessCenter(1,i) << " ";
		output << tessCenter(2,i) << " ";
		output << tessArea(i) << " ";
    }
    output.close();
}

double Cavity::compPolarizationEnergy(const std::string & potName, 
									  const std::string & chgName) 
{
	VectorXd & potVec = getFunction(potName).getVector();
	VectorXd & chgVec = getFunction(chgName).getVector();
	return potVec.dot(chgVec);
}
	
double Cavity::compPolarizationEnergy() {
	double ENN = compPolarizationEnergy("NucPot", "NucChg");
	double ENE = compPolarizationEnergy("NucPot", "EleChg");
	double EEN = compPolarizationEnergy("ElePot", "NucChg");
	double EEE = compPolarizationEnergy("ElePot", "EleChg");
//	cout << " E_ee " << EEE << " E_en " << EEN
//		 << " E_ne " << ENE << " E_nn " << ENN << endl;
	printf("E_ee = %.10E, E_en = %.10E, E_ne = %.10E, E_nn = %.10E\n", EEE, EEN, ENE, ENN);
	return 0.5 * (EEE + EEN + ENE + ENN);
}

ostream & operator<<(ostream & os, Cavity & cavity) {
	return cavity.printObject(os);
}

ostream & Cavity::printObject(ostream & os) {
	os << "Molecular cavity" << endl;
	os << "Nr. of tesserae: " << nTess;
    for(int i = 0; i < nTess; i++) {
		os << endl;
		os << i+1 << " ";
		os << tessCenter(0,i) << " ";
		os << tessCenter(1,i) << " ";
		os << tessCenter(2,i) << " ";
		os << tessArea(i);
    }
	return os;
}

/*
double Cavity::compPolarizationEnergy(std::string pot, std::string chg) {

}
*/
void Cavity::appendNewFunction(const std::string & name) {
	if(not this->isBuilt()) {
		std::cout << "Cavity not yet built!" << std::endl;
		exit(-1);
	}
	if (functions.count(name) == 0) {
		SurfaceFunction * function = new SurfaceFunction(name, nTess);
		pair<SurfaceFunctionMap::iterator, bool> retval;
		retval = functions.insert(SurfaceFunctionPair(name, function));
	}
}

void Cavity::setFunction(const std::string & name, double * values) {
	if(functions.count(name) == 0) {
		appendNewFunction(name);
	}
	SurfaceFunction * func = functions.find(name)->second;
	func->setValues(values);
}

SurfaceFunction & Cavity::getFunction(const std::string & name) {
	if(functions.count(name) == 0) {
		std::cout << "Function " << name << " does not exist" << std::endl;
		exit(-1);
	}
	SurfaceFunction * func = functions.find(name)->second;
	return * func;
}

