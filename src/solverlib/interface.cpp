#include <iostream>
#include <fstream> 
#include <string>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "PCMSolver.h"

GePolCavity *cavity;
PCMSolver *solver;

extern "C" void init_gepol_cavity_(){
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	cout << infile << endl;
	Getkw Input = Getkw(infile, false, true);
    cavity = new GePolCavity(Input);
	cavity->makeCavity();
    UniformDielectric water(78.39);
    Vacuum vacuum;
}

extern "C" void print_gepol_cavity_(){
	cout << "Cavity size" << cavity->size() << endl;
}
