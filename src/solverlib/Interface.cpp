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
#include "Solvent.h"
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

vector<Atom> Bondi = _gePolCavity->initBondi();
vector<Solvent> solventData = _solver->initSolvent();

/*

  Cavity related functions

*/

extern "C" void init_gepol_cavity_() {
	const char *infile = 0;
	enum {Explicit, Atoms, Implicit};
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	int mode = Input.getInt("Cavity<gepol>.ModeIndex");
	double scaling = Input.getDbl("Cavity<gepol>.Scaling");
	int solventIndex = Input.getInt("Medium.SolIndex");
	_gePolCavity = new GePolCavity(Input, mode, "Cavity<gepol>");
	int nNuclei;
	Matrix<double, 3, Dynamic> sphereCenter;
	VectorXd charges;
	if (Input.getStr("Cavity<gepol>.AddSpheres") == "Yes"){
		if (Input.getDbl("Cavity<gepol>.SolventRadius")) {
			_gePolCavity->setRSolv(Input.getDbl("Cavity<gepol>.SolventRadius"));
		} else {
			_gePolCavity->setRSolv(solventData[solventIndex].getSolventRadius());
		}
	} else {
		_gePolCavity->setRSolv(0.0);
	}
	init_atoms_(charges, sphereCenter);
	nNuclei = charges.size();
	_gePolCavity->setNSpheres(nNuclei);
	if ( mode == Atoms ){
		vector<int> atomsInput = Input.getIntVec("Cavity<gepol>.Atoms");
		vector<double> radiiInput = Input.getDblVec("Cavity<gepol>.Radii");
		int nsph = atomsInput.size();	
		for ( int i = 0; i < nNuclei; i++ ) {
			for ( int j = 0; j < nsph; j++ ) {
				if ( i + 1 == atomsInput[j] ) {
					charges(i) = 0;
				}
			}
		}
		for ( int i = 0; i < nNuclei; i++ ) {
			Vector3d center = sphereCenter.col(i);
			if ( charges(i) == 0 ) {
				Sphere sph(center, radiiInput[i]*scaling);
				_gePolCavity->getVectorSpheres().push_back(sph);
			} else {
				int k = charges(i) - 1;
				Sphere sph(center, Bondi[k].getAtomRadius()*scaling);
				_gePolCavity->getVectorSpheres().push_back(sph);
			}
		}
	} else if ( mode == Implicit ) {
		for ( int i = 0; i < nNuclei; i++ ) {
			int j = charges(i)-1;
			Vector3d center = sphereCenter.col(i);
			Sphere sph(center, Bondi[j].getAtomRadius()*scaling);
			_gePolCavity->getVectorSpheres().push_back(sph);
		}
	}
        _gePolCavity->makeCavity(5000, 10000000);
	_gePolCavity->initPotChg();
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
	if ( flag == 0 ) {
		sphereCenter *= ToAngstrom;
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

/*extern "C" void comp_pot_chg_pcm_(double *density, double *work, int *lwork) {
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
}*/

extern "C" void comp_pol_ene_pcm_(double * energy) {
	* energy = _cavity->compPolarizationEnergy();
}

extern "C" void print_gepol_cavity_(){
	cout << "Cavity size" << _cavity->size() << endl;
}

/*

  Solver related functions

*/

extern "C" void get_epsilon_static_(double * epsilon) {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	int solventIndex = Input.getInt("Medium.SolIndex");
        * epsilon = solventData[solventIndex].getSolventEpsStatic();
}

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
	cout << *_IEFSolver << endl;
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
extern "C" void comp_charge_(double *potential_, double *charge_, int & store) {
	enum {No, Nuclear, Electronic, Total};
	int nts = _solver->getCavitySize(); 
        switch(store) {
            case No: {
                        VectorXd potential(nts);
                        VectorXd charge(nts);
                        for (int i = 0; i < nts; i++) {
	    	                potential(i) = potential_[i];
	                }
                        charge = _solver->compCharge(potential);
	                for (int i = 0; i < nts; i++) {
                                charge_[i] = charge(i);
  	                }
                        break;
                     }
            case Nuclear: for (int i = 0; i < nts; i++) {
                              _cavity->setPot(potential_[i], Cavity::Nuclear, i); 
                          }
                          _cavity->getChg(Cavity::Nuclear) = _solver->compCharge(_cavity->getPot(Cavity::Nuclear));
                          for (int i = 0; i < nts; i++) {
                              charge_[i] = _cavity->getChg(Cavity::Nuclear, i);
                          }
                          break;
            case Electronic: for (int i = 0; i < nts; i++) {
                              _cavity->setPot(potential_[i], Cavity::Electronic, i); 
                          }
                          _cavity->getChg(Cavity::Electronic) = _solver->compCharge(_cavity->getPot(Cavity::Electronic));
                          for (int i = 0; i < nts; i++) {
                              charge_[i] = _cavity->getChg(Cavity::Electronic, i);
                          }
                          break;
            case Total: cout << "Total potential and charge coming soon!" << endl;
                        break;
        }
}
