#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "Cavity.h"
#include "Atom.h"

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

vector<Atom> Cavity::init_Bondi() {
	/*

    vector<Atom> Bondi() contains the van der Waals radii taken from
    --- A. Bondi, J. Phys. Chem. 68, 441-451 (1964) ---
  
  */

  vector<Atom> Bondi(54);
  Vector3d Origin;

  // ------------------------------------------------------------

  Origin << 0.0, 0.0, 0.0;
  
  Bondi[0] = Atom("Hydrogen", "H", 1.0, 1.20, Origin);
  Bondi[1] = Atom("Helium", "He", 2.0, 1.40, Origin);
  Bondi[2] = Atom("Lithium", "Li", 3.0, 0.0, Origin);
  Bondi[3] = Atom("Beryllium", "Be", 4.0, 0.0, Origin);
  Bondi[4] = Atom("Boron", "B", 5.0, 0.0, Origin);
  Bondi[5] = Atom("Carbon", "C", 6.0, 1.70, Origin);
  Bondi[6] = Atom("Nitrogen", "N", 7.0, 1.55, Origin);
  Bondi[7] = Atom("Oxygen", "O", 8.0, 1.52, Origin);
  Bondi[8] = Atom("Fluorine", "F", 9.0, 1.47, Origin);
  Bondi[9] = Atom("Neon", "Ne", 10.0, 1.54, Origin);
  Bondi[10] = Atom("Sodium", "Na", 11.0, 0.0, Origin);
  Bondi[11] = Atom("Magnesium", "Mg", 12.0, 0.0, Origin);
  Bondi[12] = Atom("Aluminium", "Al", 13.0, 0.0, Origin);
  Bondi[13] = Atom("Silicon", "Si", 14.0, 2.10, Origin);
  Bondi[14] = Atom("Phosphorus", "P", 15.0, 1.80, Origin);
  Bondi[15] = Atom("Sulphur", "S", 16.0, 1.80, Origin);
  Bondi[16] = Atom("Chlorine", "Cl", 17.0, 1.75, Origin);
  Bondi[17] = Atom("Argon", "Ar", 18.0, 1.88, Origin);
  Bondi[18] = Atom("Potassium", "K", 19.0, 0.0, Origin);
  Bondi[19] = Atom("Calcium", "Ca", 20.0, 0.0, Origin);
  Bondi[20] = Atom("Scandium", "Sc", 21.0, 0.0, Origin);
  Bondi[21] = Atom("Titanium", "Ti", 22.0, 0.0, Origin);
  Bondi[22] = Atom("Vanadium", "V",  23.0, 0.0, Origin);
  Bondi[23] = Atom("Chromium", "Cr", 24.0, 0.0, Origin);
  Bondi[24] = Atom("Manganese", "Mn", 25.0, 0.0, Origin);
  Bondi[25] = Atom("Iron", "Fe", 26.0, 0.0, Origin);
  Bondi[26] = Atom("Cobalt", "Co", 27.0, 0.0, Origin);
  Bondi[27] = Atom("Nickel", "Ni", 28.0, 0.0, Origin);
  Bondi[28] = Atom("Copper", "Cu", 29.0, 0.0, Origin);
  Bondi[29] = Atom("Zinc", "Zn", 30.0, 0.0, Origin);
  Bondi[30] = Atom("Gallium", "Ga", 31.0, 0.0, Origin);
  Bondi[31] = Atom("Germanium", "Ge", 32.0, 0.0, Origin);
  Bondi[32] = Atom("Arsenic", "As", 33.0, 1.85, Origin);
  Bondi[33] = Atom("Selenium", "Se", 34.0, 1.90, Origin);
  Bondi[34] = Atom("Bromine", "Br", 35.0, 1.85, Origin);
  Bondi[35] = Atom("Krypton", "Kr", 36.0, 2.02, Origin);
  Bondi[36] = Atom("Rubidium", "Rb", 37.0, 0.0, Origin);
  Bondi[37] = Atom("Strontium", "Sr", 38.0, 0.0, Origin);
  Bondi[38] = Atom("Yttrium", "Y", 39.0, 0.0, Origin);
  Bondi[39] = Atom("Zirconium", "Zr", 40.0, 0.0, Origin);
  Bondi[40] = Atom("Niobium", "Nb", 41.0, 0.0, Origin);
  Bondi[41] = Atom("Molybdenum", "Mo", 42.0, 0.0, Origin);
  Bondi[42] = Atom("Technetium", "Tc", 43.0, 0.0, Origin);
  Bondi[43] = Atom("Ruthenium", "Ru", 44.0, 0.0, Origin);
  Bondi[44] = Atom("Rhodium", "Rh", 45.0, 0.0, Origin);
  Bondi[45] = Atom("Palladium", "Pd", 46.0, 0.0, Origin);
  Bondi[46] = Atom("Silver", "Ag", 47.0, 0.0, Origin);
  Bondi[47] = Atom("Cadmium", "Cd", 48.0, 0.0, Origin);
  Bondi[48] = Atom("Indium", "In", 49.0, 0.0, Origin);
  Bondi[49] = Atom("Tin", "Sn", 50.0, 0.0, Origin);
  Bondi[50] = Atom("Antimony", "Sb", 51.0, 0.0, Origin);
  Bondi[51] = Atom("Tellurium", "Te", 52.0, 2.06, Origin);
  Bondi[52] = Atom("Iodine", "I",  53.0, 1.98, Origin);
  Bondi[53] = Atom("Xenon", "Xe", 54.0, 2.16, Origin);

  // ------------------------------------------------------------

  return Bondi;
}
