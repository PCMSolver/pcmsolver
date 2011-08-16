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
	double area = Input.getDbl("Cavity.Area");

	WaveletCavity cavity(Input);

	cavity.makeCavity();
	cout << cavity << endl;
	

}
/*
    GePolCavity cavity(Input);
	cavity.makeCavity(5000, 10000000);
	Vector3d p1(0.0, 10.1, 0.0);
    Vector3d p2(0.0, 10.2, 0.0);
    Vector3d ps(0.0,  0.0, 0.0);
    
    MetalSphere metal(20.0, 0.0, 10000.0, ps, 10.0);
    UniformDielectric water(78.39);
    Vacuum vacuum;

    PCMSolver waterSolver(vacuum, water); 

    waterSolver.buildAnisotropicMatrix(cavity);
    VectorXd potential(cavity.size());
    VectorXd charges(cavity.size());
    potential.setConstant(1.0);
    const MatrixXd &matrix = waterSolver.getPCMMatrix();
    charges = waterSolver.compCharge(potential);
	double tot = charges.sum();
	double eps = water.getEpsilon();
	cout << "Green " <<  tot << " " << tot * eps/(1.0-eps) << endl;
}
*/
