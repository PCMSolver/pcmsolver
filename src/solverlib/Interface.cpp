/*

  Interface functions implementation

*/

#include <iostream>
#include <iomanip>
#include <fstream> 
#include <string>

#include <Eigen/Dense>

#include "Constants.h"
#include "Getkw.h"
#include "taylor.hpp"
#include "SurfaceFunction.h"
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


/*

	Functions visible to host program  

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
	_solver->setSolverType(modelType);
}



extern "C" void comp_chg_pcm_(char * potName, char * chgName) {
	string potFuncName(potName);
	string chgFuncName(chgName);
	if(not _cavity->functionExists(chgFuncName)) {
		_cavity->appendNewFunction(chgFuncName);
	}
	_cavity->getFunction(chgName).clear();
	VectorXd & charge    = _cavity->getFunction(chgName).getVector();
	VectorXd & potential = _cavity->getFunction(potName).getVector();
	_solver->compCharge(potential, charge);
	double totalChg = charge.sum();
}

extern "C" void comp_pol_ene_pcm_(double * energy) {
	* energy = _cavity->compPolarizationEnergy();
}

extern "C" void get_epsilon_static_(double * epsilon) {
// This is for Gauss Theorem test on computed polarization charges
// meaningful only when there's Vacuum/UniformDielectric.
// Need to think more about this
//        * epsilon = _solver->solvent.getEpsStatic();
 	std::cout << "Not yet implemented!" << std::endl;
        exit(-1); 
}

extern "C" void get_cavity_size_(int * nts) {
	*nts = _cavity->size();
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

extern "C" void get_tess_cent_coord_(int * its, double * center) {
	Vector3d tess = _cavity->getTessCenter(*its-1);
	std::cout << tess.transpose() << std::endl;
	center[0] = tess(0);
	center[1] = tess(1);
	center[2] = tess(2);
}

extern "C" void print_gepol_cavity_(){
	cout << "Cavity size" << _cavity->size() << endl;
}

extern "C" void set_surface_function_(int * nts, double * values, char * name) {
    	int nTess = _cavity->size();
	if (nTess != *nts) {
		std::cout << "Inconsistent input" << std::endl;
	}
	std::string functionName(name);
	_cavity->setFunction(functionName, values);
}

extern "C" void get_surface_function_(int * nts, double * values, char * name) {
    	int nTess = _cavity->size();
	if (nTess != *nts) {
		std::cout << "Inconsistent input" << std::endl;
	}
	std::string functionName(name);
	VectorXd & surfaceVector = _cavity->getFunction(functionName).getVector();
	for (int i = 0; i < nTess; i++) {
		values[i] = surfaceVector(i); 
	}
}

extern "C" void add_surface_function_(char * result, double * coeff, char * part) {
	std::string resultName(result);
	std::string partName(part);
	_cavity->appendNewFunction(resultName);
	VectorXd & resultVector = _cavity->getFunction(resultName).getVector();
	VectorXd & partVector = _cavity->getFunction(partName).getVector();
	resultVector = resultVector + (*coeff) * partVector;
}

extern "C" void print_surface_function_(char * name) {
	std::string functionName(name);
	SurfaceFunction & func = _cavity->getFunction(functionName);
	std::cout << func << std::endl;
}

extern "C" bool surf_func_exists_(char * name) {
	std::string functionName(name);
	return _cavity->functionExists(functionName);
}

extern "C" void clear_surf_func_(char* name) {
	std::string functionName(name);
	SurfaceFunction & func = _cavity->getFunction(functionName);
	func.clear();
}

extern "C" void append_surf_func_(char* name) {
	std::string functionName(name);
	_cavity->appendNewFunction(functionName);
}

/*

	Functions not visible to host program

*/

void init_gepol_cavity_() {
	const char *infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	string solvent = Input.getStr("Medium.Solvent");
	Section cavSect = Input.getSect("Cavity<gepol>");
	_gePolCavity = new GePolCavity(cavSect);

	if (solvent == "Explicit") {
		_gePolCavity->setProbeRadius(Input.getDbl("Medium.ProbeRadius"));
	} else {
		SolventMap solvents = Solvent::initSolventMap();
		_gePolCavity->setProbeRadius(solvents[solvent].getRadius());
	} // No need for further checks, we have input parsing!

	VectorXd charges;
	Matrix<double, 3, Dynamic> centers;
	init_atoms_(charges, centers);
	int nAtoms = charges.size();
        _gePolCavity->setNSpheres(nAtoms);

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
}

void init_wavelet_cavity_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
    	_waveletCavity = new WaveletCavity(Input, "Cavity<wavelet>");
	_waveletCavity->makeCavity();
	_waveletCavity->readCavity("molec_dyadic.dat");
}

void init_iefsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	_IEFSolver = new IEFSolver(Medium);
	_IEFSolver->buildIsotropicMatrix(*_gePolCavity);
	cout << *_IEFSolver << endl;
}

void init_pwcsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	_PWCSolver = new PWCSolver(Medium);
	_PWCSolver->buildSystemMatrix(*_waveletCavity);
}

void init_pwlsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	_PWLSolver = new PWLSolver(Medium);
	_PWLSolver->buildSystemMatrix(*_waveletCavity);
}

void build_isotropic_matrix_() {
	_IEFSolver->buildIsotropicMatrix(*_gePolCavity);
}

void build_anisotropic_matrix_() {
	_IEFSolver->buildAnisotropicMatrix(*_gePolCavity);
}

void init_atoms_(VectorXd & charges,
	         Matrix<double, 3, Dynamic> & sphereCenter) {
	int nuclei;
	int * units;
	collect_nctot_(&nuclei);
	sphereCenter.resize(NoChange, nuclei);
	charges.resize(nuclei);
	double * chg = charges.data();
	double * centers = sphereCenter.data();
	collect_atoms_(chg, centers, units);
        if (*units == 0) {
		sphereCenter *= ToAngstrom;
        }
} 

void init_spheres_atoms_(VectorXd & charges, 
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

void init_spheres_implicit_(VectorXd & charges,	Matrix<double, 3, Dynamic> & centers)
{
	const char *infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	std::string scaling = Input.getStr("Cavity<gepol>.Scaling");
	for (int i = 0; i < charges.size(); i++) {
		vector<Atom> Bondi = Atom::initBondi();
		int index = charges(i) - 1;
		double radius = Bondi[index].getAtomRadius();
                if (scaling == "Yes") {
			radius *= Bondi[index].getAtomRadiusScaling();
                }
		Vector3d center = centers.col(i);
		Sphere sphere(center, radius);
		_gePolCavity->getSpheres().push_back(sphere);
	}
}
