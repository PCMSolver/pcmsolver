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
#include "SurfaceFunction.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "WaveletCavity.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "PCMSolver.h"
#include "IEFSolver.h"
#include "CPCMSolver.h"
#include "WEMSolver.h"
#include "PWCSolver.h"
#include "PWLSolver.h"
#include "Atom.h"
#include "Sphere.h"
#include "Solvent.h"
#include "Input.h"
#include "Interface.h"

using namespace std;
using namespace Eigen;

typedef taylor<double, 3, 1> T;

// We need globals as they must be accessible across all the functions defined
// in this interface... Could we use singletons?
Cavity        * _cavity;
GePolCavity   * _gePolCavity;
WaveletCavity * _waveletCavity;

IEFSolver * _IEFSolver;
CPCMSolver * _CPCMSolver;
PWCSolver * _PWCSolver;
PWLSolver * _PWLSolver;
PCMSolver * _solver;


/*

	Functions visible to host program  

*/

extern "C" void hello_pcm_(int * a, double * b) {
	std::cout << "Hello, PCM!" << std::endl;
	std::cout << "The integer is: " << *a << std::endl;
	std::cout << "The double is: " << *b << std::endl;
}

extern "C" void init_pcm_() {
	setupInput();
        std::string modelType = Input::TheInput().getSolverType(); 
	if (modelType == "IEFPCM") {
		_gePolCavity = initGePolCavity();
		init_iefsolver_();
		_cavity = _gePolCavity;
		_solver = _IEFSolver;
	} else if (modelType == "CPCM") {
		_gePolCavity = initGePolCavity();
		init_cpcmsolver_();
		_cavity = _gePolCavity;
		_solver = _CPCMSolver;
        } else if (modelType == "Wavelet") {
		_waveletCavity = initWaveletCavity();
		init_pwcsolver_();
		_waveletCavity->uploadPoints(_PWCSolver->getQuadratureLevel(),
									 _PWCSolver->getT_(), false);
		_cavity = _waveletCavity;
		_solver = _PWCSolver;
	} else if (modelType == "Linear") {
		_waveletCavity = initWaveletCavity();
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

// Revise this function. It's just a dirty hack now.
extern "C" void comp_pol_ene_pcm_(double * energy, int * separate_or_total) {
        if (*separate_or_total == 0) { // Using separate potentials and charges
		* energy = _cavity->compPolarizationEnergy();
        } else {
		* energy = _cavity->compPolarizationEnergy("TotPot", "TotChg") * 0.5 ;
        }
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
		Vector3d tess = _cavity->getElementCenter(i);
		centers[j] = tess(0);
		centers[j+1] = tess(1);
		centers[j+2] = tess(2);
		j += 3;
	}
	
}

extern "C" void get_tess_cent_coord_(int * its, double * center) {
	Vector3d tess = _cavity->getElementCenter(*its-1);
	std::cout << tess.transpose() << std::endl;
	center[0] = tess(0);
	center[1] = tess(1);
	center[2] = tess(2);
}

extern "C" void print_pcm_(){
        const char *infile = "@pcmsolver.inp";
        Getkw Input = Getkw(infile, false, true);
        const Section & Medium = Input.getSect("Medium");
        const string modelType = Medium.getStr("SolverType");
	if (modelType == "IEFPCM") {
		std::cout << *_IEFSolver << std::endl;
	        std::cout << *_gePolCavity << std::endl;
	} else if (modelType == "CPCM") {
		std::cout << *_CPCMSolver << std::endl;
	        std::cout << *_gePolCavity << std::endl;
        } else if (modelType == "Wavelet") {
		std::cout << *_PWCSolver << std::endl;
	        std::cout << *_waveletCavity << std::endl;
	} else if (modelType == "Linear") {
		std::cout << *_PWLSolver << std::endl;
	        std::cout << *_waveletCavity << std::endl;
	}
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

void setupInput() {
	/* Here we setup the input, meaning that we read the parsed file and store everything 
	 * it contatins inside an Input object. 
	 * This object will be unique (a Singleton) to each "run" of the module.
	 *   *** WHAT HAPPENS IN NUMERICAL OPTIMIZATIONS? ***
	 */
	Input& parsedInput = Input::TheInput();
	// The only thing we can't create immediately is the vector of spheres
	// from which the cavity is to be built.
	std::string _mode = parsedInput.getMode();
	// Get the total number of nuclei and the geometry anyway
	Eigen::VectorXd charges;
	Eigen::Matrix3Xd centers;
	initAtoms(charges, centers);

	if (_mode == "Implicit") {
		vector<Sphere> spheres = initSpheresImplicit(charges, centers);
		parsedInput.setSpheres(spheres);
	} else if (_mode == "Atoms") {
		vector<Sphere>  spheres = initSpheresAtoms(charges, centers);
		parsedInput.setSpheres(spheres);
	}
}

void initAtoms(Eigen::VectorXd & _charges, Eigen::Matrix3Xd & _sphereCenter) {
	int nuclei;
	collect_nctot_(&nuclei);
	_sphereCenter.resize(NoChange, nuclei);
	_charges.resize(nuclei);
	double * chg = _charges.data();
	double * centers = _sphereCenter.data();
	collect_atoms_(chg, centers);
} 

std::vector<Sphere> initSpheresAtoms(const Eigen::VectorXd & _charges, const Eigen::Matrix3Xd & _centers) {
	std::vector<Sphere> spheres;

	vector<int> atomsInput = Input::TheInput().getAtoms();
	vector<double> radiiInput = Input::TheInput().getRadii();
	
	for (int i = 0; i < atomsInput.size(); ++i) {
		int index = atomsInput[i] - 1; // -1 to go from human readable to machine readable
		Vector3d center = _centers.col(index);
		Sphere sph(center, radiiInput[i]);
		spheres.push_back(sph);
	}
	return spheres;
}

std::vector<Sphere> initSpheresImplicit(const Eigen::VectorXd & _charges, const Eigen::Matrix3Xd & _centers) {
	std::vector<Sphere> spheres;
	bool scaling = Input::TheInput().getScaling();

	for (int i = 0; i < _charges.size(); ++i) {
		std::vector<Atom> Bondi = Atom::initBondi();
		int index = _charges(i) - 1;
		double radius = Bondi[index].getAtomRadius();
                if (scaling) {
			radius *= Bondi[index].getAtomRadiusScaling();
                }
		Vector3d center = _centers.col(i);
		Sphere sph(center, radius);
		spheres.push_back(sph);
	}
	return spheres;
}

GePolCavity * initGePolCavity() 
{
 	double area = Input::TheInput().getArea();
	std::vector<Sphere> spheres = Input::TheInput().getSpheres();
	bool addSpheres = Input::TheInput().getAddSpheres();
	double probeRadius = Input::TheInput().getProbeRadius();
        GePolCavity * cav = new GePolCavity(area, spheres, addSpheres, probeRadius);
	return cav;
}

WaveletCavity * initWaveletCavity() {
	int patchLevel = Input::TheInput().getPatchLevel();
	std::vector<Sphere> spheres = Input::TheInput().getSpheres();
	double coarsity = Input::TheInput().getCoarsity();
	double probeRadius = Input::TheInput().getProbeRadius();
    	WaveletCavity * cav = new WaveletCavity(patchLevel, spheres, coarsity, probeRadius);
	cav->readCavity("molec_dyadic.dat");
	return cav;

}
/*
void init_wavelet_cavity_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
    	_waveletCavity = new WaveletCavity(Input, "Cavity<wavelet>");
	_waveletCavity->makeCavity();
	_waveletCavity->readCavity("molec_dyadic.dat");
}
*/

void init_iefsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	_IEFSolver = new IEFSolver(Medium);
	_IEFSolver->buildIsotropicMatrix(*_gePolCavity);
}

void init_cpcmsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	_CPCMSolver = new CPCMSolver(Medium);
	double correction = Input.getDbl("Medium<Medium>.Correction");
 	_CPCMSolver->setCorrection(correction);
	_CPCMSolver->buildIsotropicMatrix(*_gePolCavity);
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
	         Matrix3Xd & sphereCenter) {
	int nuclei;
	collect_nctot_(&nuclei);
	sphereCenter.resize(NoChange, nuclei);
	charges.resize(nuclei);
	double * chg = charges.data();
	double * centers = sphereCenter.data();
	collect_atoms_(chg, centers);
} 

void init_spheres_atoms_(VectorXd & charges, 
			Matrix3Xd & centers) {
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

void init_spheres_implicit_(VectorXd & charges,	Matrix3Xd & centers)
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
