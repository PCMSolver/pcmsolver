#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
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

void Cavity::initPotChg() {
	if (isBuilt) {
		nuclearPotential.resize(nTess);
		nuclearCharge.resize(nTess);
		electronicPotential.resize(nTess);
		electronicCharge.resize(nTess);
	} else {
		cout << "Cannot initalize Potentials and Charges" << endl;
		cout << "Cavity not yet created" << endl;
		exit(1);
	}
}

VectorXd & Cavity::getPot(const int type) {
	switch (type) 
		{
		case Nuclear :
			return nuclearPotential;
		case Electronic :
			return electronicPotential;
		default :
			cout << "Invalid request" << endl;
			exit(1);
		}
}

double Cavity::getPot(const int type, const int i) {
	switch (type) 
		{
		case Nuclear :
			return nuclearPotential(i);
		case Electronic :
			return electronicPotential(i);
		default :
			cout << "Invalid request" << endl;
			exit(1);
		}
}

VectorXd & Cavity::getChg(const int type) {
	switch (type) 
		{
		case Nuclear :
			return nuclearCharge;
		case Electronic :
			return electronicCharge;
		default :
			cout << "Invalid request" << endl;
			exit(1);
		}
}

double Cavity::getChg(const int type, const int i) {
	switch (type) 
		{
		case Nuclear :
			return nuclearCharge(i);
		case Electronic :
			return electronicCharge(i);
		default :
			cout << "Invalid request" << endl;
			exit(1);
		}
}

double Cavity::compPolarizationEnergy() {
	double EEE = electronicPotential.dot(electronicCharge);
	double EEN = electronicPotential.dot(nuclearCharge);
	double ENE = nuclearPotential.dot(electronicCharge);
	double ENN = nuclearPotential.dot(nuclearCharge);

	cout << "E_ee " << EEE << "E_en " << EEN 
		 << "E_ne " << ENE << "E_nn " << ENN << endl;
	
	return 0.5 * (EEE + EEN + ENE + ENN);
}

ostream & operator<<(ostream &os, const Cavity &cavity) {
	os << "Molecular cavity" << endl;
	os << "Nr. of tesserae: " << cavity.nTess;
    for(int i = 0; i < cavity.nTess; i++) {
		os << endl;
		os << i+1 << " ";
		os << cavity.tessCenter(0,i) << " ";
		os << cavity.tessCenter(1,i) << " ";
		os << cavity.tessCenter(2,i) << " ";
		os << cavity.tessArea(i);
    }
	return os;
}

