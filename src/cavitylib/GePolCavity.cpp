/*

GePol c++ interface and wrapper methods
written by Krzysztof Mozgawa, 2011

*/
#include <iostream>
#include <fstream> 
#include <string>
#include <vector>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "Atom.h"

GePolCavity::GePolCavity(const Getkw & Input, const string path){
	Section cavity = Input.getSect(path);
	string mode = cavity.getStr("Mode");
        averageArea = cavity.getDbl("Area");
	if ( mode == "Atoms" ){
	  //vector<int> atomsInput = cavity.getIntVec("Atoms");
	  vector<double> radiiInput = cavity.getDblVec("Radii");
	  nSpheres = radiiInput.size();
	  sphereCenter.resize(NoChange, nSpheres);
	  sphereRadius.resize(nSpheres);
	}
	else if ( mode == "Implicit" ){
	  } 
	else{
	  vector<double> spheresInput = cavity.getDblVec("Spheres");
	  nSpheres = spheresInput.size()/4; // the correctness of the size has ben checked at input parsing
	  sphereCenter.resize(NoChange, nSpheres);
	  sphereRadius.resize(nSpheres);
	  int j = 0;
	  for (int i = 0; i < nSpheres; i++) {
	    sphereCenter(0,i) = spheresInput[j];
	    sphereCenter(1,i) = spheresInput[j+1];
	    sphereCenter(2,i) = spheresInput[j+2];
	    sphereRadius(i)   = spheresInput[j+3];
	    j += 4;
	  }
	}
}


GePolCavity::GePolCavity(const Section & cavity){
	vector<double> spheresInput = cavity.getDblVec("Spheres");
	averageArea = cavity.getDbl("Area");
	nSpheres = spheresInput.size()/4; // the correctness of the size has ben checked at input parsing
    sphereCenter.resize(NoChange, nSpheres);
    sphereRadius.resize(nSpheres);
	int j = 0;
	for (int i = 0; i < nSpheres; i++) {
		sphereCenter(0,i) = spheresInput[j];
		sphereCenter(1,i) = spheresInput[j+1];
		sphereCenter(2,i) = spheresInput[j+2];
		sphereRadius(i)   = spheresInput[j+3];
		j += 4;
	}
}

bool GePolCavity::readInput(string &filename){

    ifstream input;
    input.open(filename.c_str(), fstream::in);
    if (input.eof()){
		cout << "Unexpected end of file." ;
		exit(1);
    }
    if (input.bad() != 0 ){
		cout << "Type mismatch or file corrupted." ;
		exit(1);
    }
	
    input >> nSpheres;
    sphereCenter.resize(nSpheres, NoChange);
    sphereRadius.resize(nSpheres);
    for(int i = 0; i < nSpheres; i++) {
		if (input.eof()){
			cout << "Unexpected end of file." ;
			exit(1);
		}
		if (input.bad() != 0 ){
			cout << "Type mismatch or file corrupted." ;
			exit(1);
		}
    }
    input.close();
    return false;
}

void GePolCavity::writeOutput(string &filename){

    ofstream output;
    output.open(filename.c_str(), fstream::out);

    output << nTess << endl;
    for(int i=0; i < nTess; i++) {
		output << tessCenter(0,i) << " ";
		output << tessCenter(1,i) << " ";
		output << tessCenter(2,i) << " ";
		output << tessArea(i) << " ";
		output << tessSphereCenter(0,i) << " ";
		output << tessSphereCenter(1,i) << " ";
		output << tessSphereCenter(2,i) << " ";
		output << tessRadius(i) << endl;
    }
    output.close();
}

extern"C" {
    void generatecavity_cpp_(double *xtscor, double *ytscor, double *ztscor, 
							 double *ar, double *xsphcor, double *ysphcor, 
							 double *zsphcor, double *rsph, int *nts, int *nesfp,
							 double *xe, double *ye, double *ze, double *rin, 
							 double *avgArea, double* work, int* lwork);
}

void GePolCavity::makeCavity(){
	makeCavity(10000, 10000000);
}

void GePolCavity::makeCavity(int maxts, int lwork) {

	double *xtscor  = new double[maxts];
	double *ytscor  = new double[maxts];
	double *ztscor  = new double[maxts];
	double *ar      = new double[maxts];
	double *xsphcor = new double[maxts];
	double *ysphcor = new double[maxts];
	double *zsphcor = new double[maxts];
	double *rsph    = new double[maxts];
	double *work    = new double[lwork];

	int nts;

	VectorXd xv = sphereCenter.row(0);
	VectorXd yv = sphereCenter.row(1);
	VectorXd zv = sphereCenter.row(2);

	double *xe = xv.data();
	double *ye = yv.data();
	double *ze = zv.data();

	double *rin = sphereRadius.data();
	
	generatecavity_cpp_(xtscor, ytscor, ztscor, ar, xsphcor, ysphcor, zsphcor, rsph, &nts, &nSpheres, 
						xe, ye, ze, rin, &averageArea, work, &lwork);
	
    nTess = int(nts);
    tessCenter.resize(NoChange, nTess);
    tessSphereCenter.resize(NoChange, nTess);
    tessNormal.resize(NoChange, nTess);
    tessArea.resize(nTess);
    tessRadius.resize(nTess);
    for(int i=0; i < nTess; i++){
		tessCenter(0,i) = xtscor[i];
		tessCenter(1,i) = ytscor[i];
		tessCenter(2,i) = ztscor[i];
		tessArea(i) = ar[i];
		tessSphereCenter(0,i) = xsphcor[i];
		tessSphereCenter(1,i) = ysphcor[i];
		tessSphereCenter(2,i) = zsphcor[i];
		tessRadius(i) = rsph[i];
    }

    tessNormal = tessCenter - tessSphereCenter;
    for(int i=0; i < nTess; i++){
		tessNormal.col(i) /= tessNormal.col(i).norm();
	}

	delete xtscor;
	delete ytscor;
	delete ztscor;
	delete ar;
	delete xsphcor;
	delete ysphcor;
	delete zsphcor;
	delete rsph;
	delete work;
	
	isBuilt = true;

}

vector<Atom> GePolCavity::init_Bondi() {
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

ostream & operator<<(ostream &os, const GePolCavity &cavity) {
	os << "Molecular cavity" << endl;
	os << "Nr. of spheres: " << cavity.nSpheres;
    for(int i = 0; i < cavity.nSpheres; i++) {
		os << endl;
		os << i+1 << " ";
		os << cavity.sphereCenter(0,i) << " ";
		os << cavity.sphereCenter(1,i) << " ";
		os << cavity.sphereCenter(2,i) << " ";
		os << cavity.sphereRadius(i) << " ";
    }
	return os;
}

