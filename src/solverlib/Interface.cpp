/*

  Interface functions implementation

*/

#include <iostream>
#include <iomanip>
#include <fstream> 
#include <string>

#include <Eigen/Dense>

#include "Getkw.h"
#include "taylor.hpp"
#include "Cavity.h"
#include "GePolCavity.h"
#include "WaveletCavity.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "PCMSolver.h"
#include "IEFSolver.h"
#include "WEMSolver.h"
#include "PWCSolver.h"
#include "PWLSolver.h"
#include "Atom.h"
#include "Sphere.h"
#include "Solvent.h"
#include "Interface.h"
#include "Constants.h"

using namespace std;
using namespace Eigen;

typedef taylor<double, 3, 1> T;

Cavity        * _cavity;
GePolCavity   * _gePolCavity;
WaveletCavity * _waveletCavity;

IEFSolver * _IEFSolver;
PWCSolver * _PWCSolver;
PWLSolver * _PWLSolver;
PCMSolver * _solver;

vector<Atom> Bondi = _gePolCavity->initBondi();
//vector<Solvent> solventData = _solver->initSolvent();

/*

  Cavity related functions

*/

extern "C" void init_gepol_cavity_() {
	const char *infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	string solvent = Input.getStr("Medium.Solvent");
	//	int solventIndex = Input.getInt("Medium.SolIndex");
	Section cavSect = Input.getSect("Cavity<gepol>");
	_gePolCavity = new GePolCavity(cavSect);
	//	const Keyword<double> & cavProbeRad = cavSect.getKey<double>("Cavity<gepol>.ProbeRadius");
	//	const Keyword<double> & cavityProbeRadius = Input.getKeyword<double>("Cavity<gepol>.ProbeRadius");
	//	if (cavProbeRad.isDefined()) {
	if (true) {
		std::cout << "This needs to be fixed..." << std::endl;
		_gePolCavity->setRSolv(Input.getDbl("Cavity<gepol>.ProbeRadius"));
	} else if (solvent != "Explicit") {
		SolventMap solvents = Solvent::initSolventMap();
		_gePolCavity->setRSolv(solvents[solvent]->getRadius());
	} else {
		std::cout << "Should not really get here..." << std::endl;
		exit(-1);
	}

	VectorXd charges;
	Matrix<double, 3, Dynamic> centers;
	init_atoms_(charges, centers);
	int nAtoms = charges.size();

	switch (_gePolCavity->getMode()) {
	case GePolCavity::Atoms:
		init_spheres_atoms_(charges, centers);
		break;
	case GePolCavity::Implicit:
		init_spheres_implicit_(charges, centers);
		break;
	case GePolCavity::Explicit:
		break;
	default:
		std::cout << "Case unknown" << std::endl;
		exit(-1);
	}
	_gePolCavity->makeCavity(5000, 10000000);
	//	_gePolCavity->initPotChg();
}

extern "C" void init_atoms_(VectorXd & charges,
							Matrix<double, 3, Dynamic> & sphereCenter) {
	int nuclei;
	int flag;
	collect_nctot_(&nuclei);
	sphereCenter.resize(NoChange, nuclei);
	charges.resize(nuclei);
	double * chg = charges.data();
	double * centers = sphereCenter.data();
	collect_atoms_(chg, centers, & flag);
	/*
	if ( flag == 0 ) {
		sphereCenter *= ToAngstrom;
	}
	*/
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
	nuc_pot_pcm_(&nts, _cavity->getTessCenter().data(),  
				 _cavity->getPot(Cavity::Nuclear).data());
	ele_pot_pcm_(&nts, _cavity->getTessCenter().data(),  
				 _cavity->getPot(Cavity::Electronic).data(), 
				 density, work, lwork);

	_solver->compCharge(_cavity->getPot(Cavity::Nuclear),    
					   _cavity->getChg(Cavity::Nuclear));
	_solver->compCharge(_cavity->getPot(Cavity::Electronic), 
					   _cavity->getChg(Cavity::Electronic));

	double totElChg = _cavity->getChg(Cavity::Electronic).sum();
	double totNuChg = _cavity->getChg(Cavity::Nuclear).sum();
}

extern "C" void comp_chg_pcm_() {
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

extern "C" void set_potential_(int * nts, double * potential, int * flag) {
    int nTess = _cavity->size();
	if (nTess != *nts) {
		std::cout << "Inconsistent input" << std::endl;
	}
	VectorXd & pot = _cavity->getPot(*flag);
	for (int i = 0; i < nTess; i++) {
		pot(i) = potential[i];
	}
}


/*

  Solver related functions

*/

extern "C" void init_pcm_() {
	const char *infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section & Medium = Input.getSect("Medium");
	const string modelType = Medium.getStr("SolverType");
	if (modelType == "Traditional") {
		init_gepol_cavity_();
		init_iefsolver_();
		_cavity = _gePolCavity;
		_solver = _IEFSolver;
	} else if (modelType == "Wavelet") {
		init_wavelet_cavity_();
		init_pwcsolver_();
		_waveletCavity->uploadPoints(_PWCSolver->getQuadratureLevel(),
									 _PWCSolver->getT_(), false);
		_cavity = _waveletCavity;
		_solver = _PWCSolver;
	} else if (modelType == "Linear") {
		init_wavelet_cavity_();
		init_pwlsolver_();
		_waveletCavity->uploadPoints(_PWLSolver->getQuadratureLevel(),
									 _PWLSolver->getT_(), true);
		_cavity = _waveletCavity;
		_solver = _PWLSolver;
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
	cout << *_IEFSolver << endl;
}

extern "C" void init_pwcsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	_PWCSolver = new PWCSolver(Medium);
	_PWCSolver->buildSystemMatrix(*_waveletCavity);
}

extern "C" void init_pwlsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	_PWLSolver = new PWLSolver(Medium);
	_PWLSolver->buildSystemMatrix(*_waveletCavity);
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
	_solver->compCharge(potential, charge);
	for (int i = 0; i < nts; i++) {
		charge_[i] = charge(i);
	}
}
  
extern "C" void init_spheres_atoms_(VectorXd & charges, 
						 Matrix<double, 3, Dynamic> & centers) {
	const char *infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	vector<int> atomsInput = Input.getIntVec("Cavity<gepol>.Atoms");
	vector<double> radiiInput = Input.getDblVec("Cavity<gepol>.Radii");
	for (int i = 0; i < atomsInput.size(); i++) {
		int index = atomsInput[i] - 1; // -1 to go from human readable to machine readable
		Vector3d center = centers.col(index);
		Sphere sphere(center, radiiInput[i]);
		_gePolCavity->getSpheres().push_back(sphere);
	}
}

extern "C" void init_spheres_implicit_(VectorXd & charges, 
							Matrix<double, 3, Dynamic> & centers)
{
	const char *infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	double scaling = Input.getDbl("Cavity<gepol>.Scaling");
	for (int i = 0; i < charges.size(); i++) {
		int index = charges(i) - 1;
		double radius = Bondi[index].getAtomRadius() * scaling;
		Vector3d center = centers.col(i);
		Sphere sphere(center, radius);
		_gePolCavity->getSpheres().push_back(sphere);
	}
}
