/*

  Interface functions implementation

*/

#include <iostream>
#include <iomanip>
#include <fstream> 
#include <string>

#include <Eigen/Dense>

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
#include "Atom.h"
#include "Sphere.h"
#include "Interface.h"
#include "Constants.h"

using namespace std;
using namespace Eigen;


GePolCavity * _gePolCavity;
IEFSolver * _IEFSolver;
WaveletCavity * _waveletCavity;
WEMSolver * _WEMSolver;
Cavity * _cavity;
PCMSolver * _solver;

/*

  Cavity related functions

*/

extern "C" void init_gepol_cavity_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	_gePolCavity = new GePolCavity(Input, "Cavity<gepol>");
	if ( Input.getStr("Cavity<gepol>.Mode") == "Atoms" ){
		vector<int> atomsInput = Input.getIntVec("Cavity<gepol>.Atoms");
		vector<double> radiiInput = Input.getDblVec("Cavity<gepol>.Radii");
		init_atoms_(_gePolCavity->getNSpheres(), atomsInput, 
					_gePolCavity->getSphereCenter());
	}
	else if ( Input.getStr("Cavity<gepol>.Mode") == "Implicit" ){
		init_implicit_(_gePolCavity->getSphereRadius(), _gePolCavity->getSphereCenter());
	}
	_gePolCavity->makeCavity(5000, 10000000);
	_gePolCavity->initPotChg();
}

extern "C" void init_atoms_(int nSpheres, vector<int> & atomsInput, 
							Matrix<double, 3, Dynamic> & sphereCenter){
	double * centers = sphereCenter.data();
	int * idx = new int[nSpheres];
	int flag;
	for ( int i = 0; i < nSpheres; i++ ) {
		idx[i] = atomsInput[i];
	}
	collect_atoms_(& nSpheres, idx, centers, & flag);
	delete idx;
	if ( flag == 0 ) {
		sphereCenter *= ToAngstrom;
	}
}

extern "C" void init_implicit_(VectorXd & sphereRadius, 
							   Matrix<double, 3, Dynamic> & sphereCenter) {
	vector<Atom> Bondi = _gePolCavity->init_Bondi();
	VectorXd charges;
	int nuclei;
	int flag;
	collect_nctot_(&nuclei);
	_gePolCavity->setNSpheres(nuclei);
	sphereRadius.resize(nuclei);
	sphereCenter.resize(NoChange, nuclei);
	charges.resize(nuclei);
	double * chg = charges.data();
	double * centers = sphereCenter.data();
	collect_implicit_(chg, centers, & flag);
	if ( flag == 0 ) {
		sphereCenter *= ToAngstrom;
	}
	for ( int i = 0; i < nuclei; i++ ) {
		int j = charges(i)-1;
		sphereRadius(i) = Bondi[j].getAtomRadius();
	}
} 

extern "C" void init_wavelet_cavity_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
    _waveletCavity = new WaveletCavity(Input, "Cavity<wavelet>");
	_waveletCavity->makeCavity();
	_waveletCavity->readCavity("molec_dyadic.dat");
}

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

extern "C" void print_gepol_cavity_(){
	cout << "Cavity size" << _cavity->size() << endl;
}

/*

  Solver related functions

*/


extern "C" void init_pcm_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const string modelType = Input.getStr("SolverType");
	if (modelType == "Traditional") {
		init_gepol_cavity_();
		init_iefsolver_();
		_cavity = _gePolCavity;
		_solver = _IEFSolver;
	} else if (modelType == "Wavelet") {
		init_wavelet_cavity_();
		init_wemsolver_();
		_waveletCavity->uploadPoints(_WEMSolver->getQuadratureLevel(),
									 _WEMSolver->getT_());
		_cavity = _waveletCavity;
		_solver = _WEMSolver;
	} 
	_cavity->initPotChg();
	_solver->setSolverType(modelType);
}

/*
extern "C" void init_gepol_cavity_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
    _gePolCavity = new GePolCavity(Input, "Cavity<gepol>");
	_gePolCavity->makeCavity(5000, 10000000);
}
*/



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
  
