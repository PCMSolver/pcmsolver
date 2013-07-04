/*

  Interface functions implementation

*/

#include <iostream>
#include <iomanip>
#include <fstream> 
#include <string>
#include <stdexcept>

#include <Eigen/Dense>

#include "Getkw.h"
#include "SurfaceFunction.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "WaveletCavity.h"
#include "CavityFactory.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "GreensFunctionFactory.h"
#include "PCMSolver.h"
#include "IEFSolver.h"
#include "CPCMSolver.h"
#include "WEMSolver.h"
#include "PWCSolver.h"
#include "PWLSolver.h"
#include "SolverFactory.h"
#include "Atom.h"
#include "Sphere.h"
#include "Solvent.h"
#include "Input.h"
#include "Interface.h"

using namespace std;
using namespace Eigen;

typedef std::map<std::string, SurfaceFunction *> SurfaceFunctionMap;

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
//        std::string modelType = Input::TheInput().getSolverType();
        initCavity();
	initSolver();
	/*
	if (modelType == "IEFPCM") {
//		_gePolCavity = initGePolCavity();
		initSolver();
//		init_iefsolver_();
//		_cavity = _gePolCavity;
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
	_solver->setSolverType(modelType);*/
}



extern "C" void comp_chg_pcm_(char * potName, char * chgName) 
{
	string potFuncName(potName);
	string chgFuncName(chgName);

	// Get the SurfaceFunctionMap
	SurfaceFunctionMap& functions = SurfaceFunction::TheMap();
	
	// Get the proper iterators
	SurfaceFunctionMap::const_iterator iter_pot = functions.find(potFuncName);
	SurfaceFunctionMap::const_iterator iter_chg = functions.find(chgFuncName);

	if( iter_chg == functions.end() ) 
	{
		SurfaceFunction * func = new SurfaceFunction(chgFuncName, _cavity->size());
		// The SurfaceFunction automagically registers itself in the map.
		// We must update iter_chg, otherwise it will point somewhere else
		++iter_chg;
	} 
	// If it already exists there's no problem, we will pass a reference to its values to
	// _solver->compCharge(*, *) so they will be automagically updated!

	VectorXd & charge = iter_chg->second->getVector();
	VectorXd & potential = iter_pot->second->getVector();
	
	_solver->compCharge(potential, charge);
	double totalChg = charge.sum();
}

// Revise this function. It's just a dirty hack now.
extern "C" void comp_pol_ene_pcm_(double * energy, int * separate_or_total) 
{
	SurfaceFunctionMap& functions = SurfaceFunction::TheMap();
        if (*separate_or_total == 0) 
	{ // Using separate potentials and charges
		SurfaceFunctionMap::const_iterator iter_nuc_pot = functions.find("NucPot");
		SurfaceFunctionMap::const_iterator iter_nuc_chg = functions.find("NucChg");
		SurfaceFunctionMap::const_iterator iter_ele_pot = functions.find("ElePot");
		SurfaceFunctionMap::const_iterator iter_ele_chg = functions.find("EleChg");

		double UNN = (*iter_nuc_pot->second) *  (*iter_nuc_chg->second);
		double UEN = (*iter_ele_pot->second) *  (*iter_nuc_chg->second);
		double UNE = (*iter_nuc_pot->second) *  (*iter_ele_chg->second);
		double UEE = (*iter_ele_pot->second) *  (*iter_ele_chg->second);
		
		printf("U_ee = %.10E, U_en = %.10E, U_ne = %.10E, U_nn = %.10E\n", UEE, UEN, UNE, UNN);

		*energy = 0.5 * ( UNN + UEN + UNE + UEE );
        } 
	else 
	{
		SurfaceFunctionMap::const_iterator iter_pot = functions.find("TotPot");
		SurfaceFunctionMap::const_iterator iter_chg = functions.find("TotChg");
	
		*energy = (*iter_pot->second) * (*iter_chg->second) * 0.5;
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
//	        std::cout << *_gePolCavity << std::endl;
		std::cout << *_cavity << std::endl;
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

extern "C" void set_surface_function_(int * nts, double * values, char * name)
{
	int nTess = _cavity->size();
	if ( nTess != *nts )
		throw std::runtime_error("You are trying to allocate a SurfaceFunction bigger than the cavity!");

	std::string functionName(name);
	// Here we check whether the function exists already or not
	SurfaceFunctionMap & functions = SurfaceFunction::TheMap();
	SurfaceFunctionMap::const_iterator iter = functions.find(functionName);
	if ( iter == functions.end() )
	{	// If not create it!
	        SurfaceFunction * func = new SurfaceFunction(functionName, *nts, values);
		// The SurfaceFunction automagically registers itself in the map.
	}
	else
	{	// If yes just update the values!
	        iter->second->setValues(values);
	}
}

extern "C" void get_surface_function_(int * nts, double * values, char * name) 
{
    	int nTess = _cavity->size();
	if ( nTess != *nts ) 
		throw std::runtime_error("You are trying to access a SurfaceFunction bigger than the cavity!");
	
	std::string functionName(name);
	
	SurfaceFunctionMap & functions = SurfaceFunction::TheMap();
	SurfaceFunctionMap::const_iterator iter = functions.find(functionName);
	
	Eigen::VectorXd surfaceVector = iter->second->getVector();
	
	for ( int i = 0; i < nTess; ++i ) 
		values[i] = surfaceVector(i); 
}

extern "C" void add_surface_function_(char * result, double * coeff, char * part) 
{
	std::string resultName(result);
	std::string partName(part);

	append_surf_func_(result);
	
	SurfaceFunctionMap & functions = SurfaceFunction::TheMap();

	SurfaceFunctionMap::const_iterator iter_part = functions.find(partName);
	SurfaceFunctionMap::const_iterator iter_result = functions.find(resultName);

	// Using iterators and operator overloading: so neat!
	(*iter_result->second) += (*coeff) * (*iter_part->second);
}

extern "C" void print_surface_function_(char * name) 
{
	std::string functionName(name);

	SurfaceFunctionMap & functions = SurfaceFunction::TheMap();
	SurfaceFunctionMap::const_iterator iter = functions.find(name);

	std::cout << *(iter->second) << std::endl;
}

extern "C" bool surf_func_exists_(char * name) 
{
	std::string functionName(name);

	SurfaceFunctionMap & functions = SurfaceFunction::TheMap();
	SurfaceFunctionMap::const_iterator iter = functions.find(name);

	return iter != functions.end();
}

extern "C" void clear_surf_func_(char* name) 
{
	std::string functionName(name);

	SurfaceFunctionMap & functions = SurfaceFunction::TheMap();
	SurfaceFunctionMap::const_iterator iter = functions.find(name);

	iter->second->clear();
}

extern "C" void append_surf_func_(char* name) 
{
	int nTess = _cavity->size();
	std::string functionName(name);

	SurfaceFunctionMap & functions = SurfaceFunction::TheMap();

	SurfaceFunctionMap::const_iterator iter = functions.find(functionName);
	// Check if function already exists in the map
	if ( iter == functions.end() )
	{	// If not create it
		SurfaceFunction * func = new SurfaceFunction(functionName, nTess);
		// The SurfaceFunction automagically registers itself in the map.
	} // What happens if it is already in the map?
}

/*

	Functions not visible to host program

*/

void setupInput() {
	/* Here we setup the input, meaning that we read the parsed file and store everything 
	 * it contatins inside an Input object. 
	 * This object will be unique (a Singleton) to each "run" of the module.
	 *   *** WHAT HAPPENS IN NUMERICAL GEOMETRY OPTIMIZATIONS? ***
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

void initCavity()
{
	// Get the input data for generating the cavity
	std::string cavityType = Input::TheInput().getCavityType();
 	double area = Input::TheInput().getArea();
	std::vector<Sphere> spheres = Input::TheInput().getSpheres();
	bool addSpheres = Input::TheInput().getAddSpheres();
	double probeRadius = Input::TheInput().getProbeRadius();
	int patchLevel = Input::TheInput().getPatchLevel();
	double coarsity = Input::TheInput().getCoarsity();

	// Get the right cavity from the Factory
	_cavity = CavityFactory::TheCavityFactory().createCavity(cavityType, spheres, area, probeRadius, addSpheres, patchLevel, coarsity);
	// Our use of inheritance breaks Liskov Principle, the Factory is pretty useless... 
	// UNLESS one uses RTTI (which we do) although I don't know if it's a sign of ugly design strategy... 
}

void initSolver()
{
	GreensFunctionFactory & factory = GreensFunctionFactory::TheGreensFunctionFactory();
	// Get the input data for generating the inside & outside Green's functions
	// INSIDE
	double epsilon = Input::TheInput().getEpsilonInside();
	std::string greenType = Input::TheInput().getGreenInsideType();
	int greenDer = Input::TheInput().getDerivativeInsideType();

	GreensFunction * gfInside = factory.createGreensFunction(greenType, greenDer);
	
	// OUTSIDE, reuse the variables holding the parameters for the Green's function inside.
	epsilon = Input::TheInput().getEpsilonOutside();
	greenType = Input::TheInput().getGreenOutsideType();
	greenDer = Input::TheInput().getDerivativeOutsideType();
	
	GreensFunction * gfOutside = factory.createGreensFunction(greenType, greenDer, epsilon);
	// And all this to finally create the solver! 
	std::string modelType = Input::TheInput().getSolverType();
	double correction = Input::TheInput().getCorrection();
	int eqType = Input::TheInput().getEquationType();
	_solver = SolverFactory::TheSolverFactory().createSolver(modelType, gfInside, gfOutside, correction, eqType);
	_solver->buildSystemMatrix(*_cavity);
/*	if (modelType == "IEFPCM") {
		_IEFSolver = new IEFSolver(gfInside, gfOutside);
		std::cout << "solver done" << std::endl;
		_IEFSolver->buildSystemMatrix(*_cavity);
		std::cout << "system matrix done" << std::endl;
		_solver = _IEFSolver;
	} else if (modelType == "CPCM") {
		_IEFSolver = new IEFSolver(gfInside, gfOutside);
		double correction = Input::TheInput().getCorrection();
        	_CPCMSolver->setCorrection(correction);
		_CPCMSolver->buildSystemMatrix(*_cavity);
		_solver = _CPCMSolver;
        } else*/ if (modelType == "Wavelet") {
		_waveletCavity = initWaveletCavity();
		_PWCSolver = new PWCSolver(gfInside, gfOutside);
		_PWCSolver->buildSystemMatrix(*_waveletCavity);
		_waveletCavity->uploadPoints(_PWCSolver->getQuadratureLevel(), _PWCSolver->getT_(), false); // WTF is happening here???
		_cavity = _waveletCavity;
		_solver = _PWCSolver;
	} else if (modelType == "Linear") {
		_waveletCavity = initWaveletCavity();
		_PWLSolver = new PWLSolver(gfInside, gfOutside);
		_PWLSolver->buildSystemMatrix(*_waveletCavity);
		_waveletCavity->uploadPoints(_PWLSolver->getQuadratureLevel(),_PWLSolver->getT_(), true); // WTF is happening here???
		_cavity = _waveletCavity;
		_solver = _PWLSolver;
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
        GePolCavity * cav = new GePolCavity(spheres, area, probeRadius, addSpheres);
	return cav;
}

WaveletCavity * initWaveletCavity() {
	int patchLevel = Input::TheInput().getPatchLevel();
	std::vector<Sphere> spheres = Input::TheInput().getSpheres();
	double coarsity = Input::TheInput().getCoarsity();
	double probeRadius = Input::TheInput().getProbeRadius();
    	WaveletCavity * cav = new WaveletCavity(spheres, probeRadius, patchLevel, coarsity);
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
//	_IEFSolver = new IEFSolver(Medium);
//	_IEFSolver->buildIsotropicMatrix(*_gePolCavity);
//	_IEFSolver->buildSystemMatrix(*_cavity);
}

void init_cpcmsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
//	_CPCMSolver = new CPCMSolver(Medium);
//	double correction = Input.getDbl("Medium<Medium>.Correction");
// 	_CPCMSolver->setCorrection(correction);
//	_CPCMSolver->buildIsotropicMatrix(*_gePolCavity);
}

void init_pwcsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	//_PWCSolver = new PWCSolver(Medium);
	//_PWCSolver->buildSystemMatrix(*_waveletCavity);
}

void init_pwlsolver_() {
	const char *infile = 0;
	infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
	const Section &Medium = Input.getSect("Medium<Medium>");
	//_PWLSolver = new PWLSolver(Medium);
	//_PWLSolver->buildSystemMatrix(*_waveletCavity);
}

void build_isotropic_matrix_() {
//	_IEFSolver->buildIsotropicMatrix(*_gePolCavity);
}

void build_anisotropic_matrix_() {
//	_IEFSolver->buildAnisotropicMatrix(*_gePolCavity);
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
