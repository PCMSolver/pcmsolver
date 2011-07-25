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

//      Subroutine PotExpVal(Density, Centers, Nts, Potential, Work, 
//     $                     LWork)
extern "C" void electron_pot_(double *density, double* centers, int *nts, 
							 double *potential, double *work, double *lwork)

extern "C" void nuclear_pot_(double* centers, int *nts, double *potential)

//      Subroutine Fock_PCMModule(Fock, Centers, Nts, Charges, Work, 
//     $                     LWork)
extern "C" void fock_pcm_module_(double *fock, double* centers, int *nts, 
								 double *charges, double *work, double *lwork);

extern "C" void pcm_scf_(double *fock, double *energy, double *density, 
						 double *work, int *lwork) {
	VectorXd ElectronPotential(cavity->size());
	VectorXd NuclearPotential(cavity->size());

	nuclear_pot_(cavity->getCenterTess.data(), cavity->size(), NuclearPotential.data())
	electron_pot_(density, cavity->getCenterTess.data(), cavity->size(), 
				 ElectronPotential.data(), work, lwork);
	VectorXd NuclearCharge = solver.compCharge(NuclearPotential);
	VectorXd ElectronCharge = solver.compCharge(ElectronPotential);
	VectorXd TotalCharge = ElectronCharge + NuclearCharge;
	Eee = ElectronCharge * ElectronPotential;
	Een = ElectronCharge * NuclearPotential;
	Ene = NuclearCharge  * ElectronPotential;
	Enn = NuclearCharge  * NuclearPotential;
	Etot = 0.5 * (Eee + Een + Ene + Enn);
	fock_pcm_module_(fock, cavity->getCenterTess.data(), cavity->size(),  
					 ElectronCharge.data(), work, lwork);
}


extern "C" void init_gepol_cavity_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
    cavity = new GePolCavity(Input);
	cavity->makeCavity(5000, 10000000);
}

extern "C" void init_pcmsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	cout << Medium << endl;
	solver = new PCMSolver(Medium);
}

extern "C" void build_anisotropic_matrix_() {
	solver->buildAnisotropicMatrix(*cavity);
}

//copying mechanism of the following routine needs to be revised
extern "C" void comp_charge_(double *potential_, double *charge_) {
	int nts = solver->getCavitySize(); 
	VectorXd potential(nts);
	VectorXd charge(nts);
	for (int i = 0; i < nts; i++) {
		potential(i) = potential_[i];
	}
	charge = solver->compCharge(potential);
	for (int i = 0; i < nts; i++) {
		charge_[i] = charge(i);
	}
}

extern "C" void print_gepol_cavity_(){
	cout << "Cavity size" << cavity->size() << endl;
}

