#include <iostream>
#include <iomanip>
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
extern "C" void ele_pot_pcm_(double *density, double* centers, int *nts, 
							  double *potential, double *work, int *lwork);

extern "C" void nuc_pot_pcm_(double* centers, int *nts, double *potential);

//      Subroutine Fock_PCMModule(Fock, Centers, Nts, Charges, Work, 
//     $                     LWork)
extern "C" void fock_pcm_(double *fock, double* centers, int *nts, 
								 double *charges, double *work, int *lwork);

extern "C" void get_cavity_size_(int * nts) {
	*nts = cavity->size();
}

extern "C" void get_surface_charge_(double * charge) {
	for (int i = 0; i < cavity->size(); i++) {
		charge[i] = cavity->getChg(Cavity::Nuclear, i) + 
			        cavity->getChg(Cavity::Electronic, i);
	}
}

extern "C" void get_tess_centers_(double * centers) {
	int j = 0;
	for (int i = 0; i < cavity->size(); i++) {
		Vector3d tess = cavity->getTessCenter(i);
		centers[j] = tess(0);
		centers[j+1] = tess(1);
		centers[j+2] = tess(2);
		j += 3;
	}
}

extern "C" void comp_pot_chg_pcm_(double *density, double *work, int *lwork) {
	int nts = cavity->size();
	nuc_pot_pcm_(cavity->getTessCenter().data(), &nts, 
				 cavity->getPot(Cavity::Nuclear).data());
	ele_pot_pcm_(density, cavity->getTessCenter().data(), &nts, 
				 cavity->getPot(Cavity::Electronic).data(), work, lwork);

	solver->compCharge(cavity->getPot(Cavity::Nuclear),    
					   cavity->getChg(Cavity::Nuclear));
	solver->compCharge(cavity->getPot(Cavity::Electronic), 
					   cavity->getChg(Cavity::Electronic));

	double totElChg = cavity->getChg(Cavity::Electronic).sum();
	double totNuChg = cavity->getChg(Cavity::Nuclear).sum();

	cout << "total charges" << " " << totElChg << " " << totNuChg << endl;
}

extern "C" void comp_pol_ene_pcm_(double * energy) {
	* energy = cavity->compPolarizationEnergy();
}

extern "C" void init_gepol_cavity_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
    cavity = new GePolCavity(Input);
	cout << "before make cavity" << endl;
	cavity->makeCavity(5000, 10000000);
	cout << "after make cavity" << endl;
	cavity->initPotChg();
}

extern "C" void init_pcmsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	solver = new PCMSolver(Medium);
	solver->buildIsotropicMatrix(*cavity);
}

extern "C" void build_anisotropic_matrix_() {
	solver->buildAnisotropicMatrix(*cavity);
}

extern "C" void build_isotropic_matrix_() {
	solver->buildIsotropicMatrix(*cavity);
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

