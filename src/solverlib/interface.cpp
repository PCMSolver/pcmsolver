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
#include "WaveletCavity.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "PCMSolver.h"
#include "IEFSolver.h"
#include "WEMSolver.h"

GePolCavity * _gePolCavity;
IEFSolver * _IEFSolver;
WaveletCavity * _waveletCavity;
WEMSolver * _WEMSolver;
Cavity * _cavity;
PCMSolver * _solver;

//      Subroutine PotExpVal(Density, Centers, Nts, Potential, Work, 
//     $                     LWork)
extern "C" void ele_pot_pcm_(double *density, double* centers, int *nts, 
							  double *potential, double *work, int *lwork);

extern "C" void nuc_pot_pcm_(double* centers, int *nts, double *potential);

//      Subroutine Fock_PCMModule(Fock, Centers, Nts, Charges, Work, 
//     $                     LWork)
extern "C" void fock_pcm_(double *fock, double* centers, int *nts, 
								 double *charges, double *work, int *lwork);

extern "C" void init_gepol_cavity_();

extern "C" void init_iefsolver_();

extern "C" void init_wavelet_cavity_();

extern "C" void init_wemsolver_();

extern "C" void init_pcm_();

extern "C" void get_cavity_size_(int * nts) {
	*nts = _cavity->size();
}

extern "C" void get_total_surface_charge_(double * charge) {
	for (int i = 0; i < _cavity->size(); i++) {
		charge[i] = _cavity->getChg(Cavity::Nuclear, i) + 
			        _cavity->getChg(Cavity::Electronic, i);
	}
}

extern "C" void get_nuclear_surface_charge_(double * charge) {
	for (int i = 0; i < _cavity->size(); i++) {
		charge[i] = _cavity->getChg(Cavity::Nuclear, i);
	}
}

extern "C" void get_electronic_surface_charge_(double * charge) {
	for (int i = 0; i < _cavity->size(); i++) {
		charge[i] = _cavity->getChg(Cavity::Electronic, i);
	}
}

extern "C" void get_tess_centers_(double * centers) {
	int j = 0;
	for (int i = 0; i < _cavity->size(); i++) {
		Vector3d tess = _cavity->getTessCenter(i);
		centers[j] = tess(0);
		centers[j+1] = tess(1);
		centers[j+2] = tess(2);
		j += 3;
	}
	
}

extern "C" void comp_pot_chg_pcm_(double *density, double *work, int *lwork) {
	int nts = _cavity->size();
	nuc_pot_pcm_(_cavity->getTessCenter().data(), &nts, 
				 _cavity->getPot(Cavity::Nuclear).data());
	ele_pot_pcm_(density, _cavity->getTessCenter().data(), &nts, 
				 _cavity->getPot(Cavity::Electronic).data(), work, lwork);

	_solver->compCharge(_cavity->getPot(Cavity::Nuclear),    
					   _cavity->getChg(Cavity::Nuclear));
	_solver->compCharge(_cavity->getPot(Cavity::Electronic), 
					   _cavity->getChg(Cavity::Electronic));

	double totElChg = _cavity->getChg(Cavity::Electronic).sum();
	double totNuChg = _cavity->getChg(Cavity::Nuclear).sum();
}

extern "C" void comp_pol_ene_pcm_(double * energy) {
	* energy = _cavity->compPolarizationEnergy();
}

extern "C" void init_pcm_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const string modelType = Input.getStr("SolverType");
	if (modelType == "Traditional") {
		cout << "Traditional chosen" << endl;
		init_gepol_cavity_();
		cout << "Cavity made" << endl;
		init_iefsolver_();
		cout << "Solver made" << endl;
		_cavity = _gePolCavity;
		cout << "cavity pointer" << endl;
		_solver = _IEFSolver;
		cout << "solver pointer" << endl;
	} else if (modelType == "Wavelet") {
		init_wavelet_cavity_();
		init_wemsolver_();
		_waveletCavity->uploadPoints(_WEMSolver->getQuadratureLevel(),
									 _WEMSolver->getT_());
		_cavity = _waveletCavity;
		_solver = _WEMSolver;
	} 
	_solver->setSolverType(modelType);
}

extern "C" void init_gepol_cavity_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
    _gePolCavity = new GePolCavity(Input, "Cavity<gepol>");
	_gePolCavity->makeCavity(5000, 10000000);
	_gePolCavity->initPotChg();
}

extern "C" void init_wavelet_cavity_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
    _waveletCavity = new WaveletCavity(Input, "Cavity<wavelet>");
	_waveletCavity->makeCavity();
	_waveletCavity->readCavity("molec_dyadic.dat");
	_waveletCavity->initPotChg();
}

extern "C" void init_iefsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	_IEFSolver = new IEFSolver(Medium);
	_IEFSolver->buildIsotropicMatrix(*_gePolCavity);
}

extern "C" void init_wemsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	_WEMSolver = new WEMSolver(Medium);
	_WEMSolver->uploadCavity(*_waveletCavity);
	_WEMSolver->constructSystemMatrix();
}

extern "C" void build_anisotropic_matrix_() {
	_IEFSolver->buildAnisotropicMatrix(*_gePolCavity);
}

extern "C" void build_isotropic_matrix_() {
	_IEFSolver->buildIsotropicMatrix(*_gePolCavity);
}

//copying mechanism of the following routine needs to be revised
extern "C" void comp_charge_(double *potential_, double *charge_) {
	int nts = _solver->getCavitySize(); 
	VectorXd potential(nts);
	VectorXd charge(nts);
	for (int i = 0; i < nts; i++) {
		potential(i) = potential_[i];
	}
	charge = _solver->compCharge(potential);
	for (int i = 0; i < nts; i++) {
		charge_[i] = charge(i);
	}
}

extern "C" void print_gepol_cavity_(){
	cout << "Cavity size" << _cavity->size() << endl;
}

