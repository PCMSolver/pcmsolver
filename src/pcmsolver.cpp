#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

extern "C"{
#include "vector3.h"
#include "sparse2.h"
#include "intvector.h"
#include "basis.h"
#include "WEM.h"
#include "read_points.h"
#include "vector2.h"
#include "interpolate.h"
#include "topology.h"
#include "kern.h"
#include "compression.h"
#include "postproc.h"
#include "WEMRHS.h"
#include "WEMPCG.h"
#include "WEMPGMRES.h"
#include "dwt.h"
#include "cubature.h"
#include "gauss_square.h"
#include "constants.h"
}



#include "Getkw.h"
#include "taylor.hpp"
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
#include "PWCSolver.h"
#include "PWLSolver.h"

void IEFTest(const Getkw & input);
template <class SolverClass> void WEMTest(const Getkw & input);

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
	Getkw input = Getkw(infile, false, true);
	std::string solverType = input.getStr("Medium.SolverType");
	if (solverType == "IEFPCM") {
		IEFTest(input);
	} else if (solverType == "Wavelet") {
		std::cout << "PWC solver...." << std::endl;
		WEMTest<PWCSolver>(input);
	} else if (solverType == "Linear") {
		std::cout << "PWL solver...." << std::endl;
		WEMTest<PWLSolver>(input);
	}
}

void IEFTest(const Getkw & input) {
	std::cout << "IEFTest NYI" << std::endl;
}

template<class SolverClass>
void WEMTest(const Getkw & input) {
	const Section & Medium = input.getSect("Medium");
	const Section & WaveletCavitySection = input.getSect("Cavity<wavelet>");
    WaveletCavity wavcav(WaveletCavitySection);
	wavcav.makeCavity();
	string wavcavFile = "molec_dyadic.dat";
	wavcav.readCavity(wavcavFile);
    SolverClass waveletSolver(Medium);
	waveletSolver.buildSystemMatrix(wavcav);
   	wavcav.uploadPoints(waveletSolver.getQuadratureLevel(), 
						waveletSolver.getT_(), waveletSolver.isPWL());
	wavcav.compFakePotential();
	waveletSolver.buildCharge(wavcav, "NucPot", "NucChg");
	double energy = wavcav.compPolarizationEnergy("NucPot","NucChg");
	cout << "Energy: " << energy << endl;
}

template void WEMTest<PWCSolver>(const Getkw & input);
template void WEMTest<PWLSolver>(const Getkw & input);

/*
	double factor = 78.39/77.39;
	double q = factor * (wavcav.getChg(Cavity::Nuclear)).sum();
	cout << " charges computed " << q << endl;
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

