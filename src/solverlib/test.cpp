#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "MetalSphere.h"
#include "GreensFunctionSum.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "WaveletCavity.h"
#include "PCMSolver.h"
#include "IEFSolver.h"
#include "WEMSolver.h"

int main(int argc, char** argv){

	const char *infile = 0;
	if (argc == 1) {
		infile = "STDIN";
	} else if (argc == 2) {
		infile = argv[1];
	} else {
		cout << "Invalid nr. of arguments" << endl;
		exit(1);
	}
	
	Getkw Input = Getkw(infile, false, true);

	int printl = Input.getInt("PRINTL");

	//	WaveletCavity wavcav(Input);

	//	cavity.makeCavity();
	//	cout << cavity << endl;
	
	Vector3d p1(0.0, 10.1, 0.0);
    Vector3d p2(0.0, 10.2, 0.0);
    Vector3d ps(0.0,  0.0, 0.0);
    
	const Section &Medium = Input.getSect("Medium");
	const Section &WaveletCavitySection = Input.getSect("Cavity<wavelet>");
	const Section &GepolCavitySection = Input.getSect("Cavity<gepol>");

    GePolCavity cavity(GepolCavitySection);
    WaveletCavity wavcav(WaveletCavitySection);

	cavity.makeCavity();
	wavcav.makeCavity();

	string wavcavFile = "molec_dyadic.dat";

	wavcav.readCavity(wavcavFile);

	cout << wavcav << endl;

    WEMSolver waveletSolver(Medium);
	cout << "wavelet solver initialized" << endl;
	waveletSolver.uploadCavity(wavcav);
	waveletSolver.constructSystemMatrix();
	cout << "system matix built" << endl;

	/*
    IEFSolver waterSolver(Medium); 
    waterSolver.buildAnisotropicMatrix(cavity);
    VectorXd potential(cavity.size());
    VectorXd charges(cavity.size());
    potential.setConstant(1.0);
    const MatrixXd &matrix = waterSolver.getPCMMatrix();
    charges = waterSolver.compCharge(potential);
	double tot = charges.sum();
	cout << "Cavity " << endl << cavity << endl;
	cout << "Cavity " << endl << cavity.size() << endl;
	cout << "Green " <<  tot << " " << tot * 78.39/(1.0-78.39) << endl;
	*/
}

