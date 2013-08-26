/*

  Interface functions implementation

*/
#include "Interface.hpp"

#include <iostream>
#include <map>
#include <string>
#include <stdexcept>
#include <utility>

#include "Config.hpp"

#include <Eigen/Dense>

// Core classes
//    1. Cavities
#include "Cavity.hpp"
#include "GePolCavity.hpp"
#include "WaveletCavity.hpp"
//    2. Green's functions
#include "GreensFunction.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
//    3. Solvers
#include "IEFSolver.hpp"
#include "CPCMSolver.hpp"
#include "PCMSolver.hpp"
#include "PWCSolver.hpp"
#include "PWLSolver.hpp"
#include "WEMSolver.hpp"
// The factories
#include "CavityFactory.hpp"
#include "GreensFunctionFactory.hpp"
#include "SolverFactory.hpp"
// Helper classes
#include "Atom.hpp"
#include "Input.hpp"
#include "Solvent.hpp"
#include "Sphere.hpp"
#include "SurfaceFunction.hpp"

typedef std::map<std::string, SurfaceFunction *> SurfaceFunctionMap;
typedef std::pair<std::string, SurfaceFunction *> SurfaceFunctionPair;

// We need globals as they must be accessible across all the functions defined in this interface...
// The final objective is to have only a pointer to Cavity and a pointer to PCMSolver (our abstractions)
// then maybe manage them through "objectification" of this interface.
Cavity        * _cavity;
WaveletCavity * _waveletCavity;

PWCSolver * _PWCSolver;
PWLSolver * _PWLSolver;
PCMSolver * _solver;

SurfaceFunctionMap functions;

/*

	Functions visible to host program  

*/

extern "C" void hello_pcm_(int * a, double * b) 
{
	std::cout << "Hello, PCM!" << std::endl;
	std::cout << "The integer is: " << *a << std::endl;
	std::cout << "The double is: " << *b << std::endl;
}

extern "C" void init_pcm_() 
{
	setupInput();
        initCavity();
	initSolver();
}

extern "C" void tear_down_pcm_()
{// Delete all the global pointers, maybe in a more refined way...

	for ( SurfaceFunctionMap::iterator iter = functions.begin(); iter != functions.end(); ++iter )
	{// Delete all SurfaceFunction
		delete (iter->second);
		functions.erase(iter);
	}
	functions.clear();

	delete _cavity;
	delete _waveletCavity;
	delete _PWCSolver;
	delete _PWLSolver;
	delete _solver;
}

extern "C" void comp_chg_pcm_(char * potName, char * chgName) 
{
	std::string potFuncName(potName);
	std::string chgFuncName(chgName);

	// Get the proper iterators
	SurfaceFunctionMap::const_iterator iter_pot = functions.find(potFuncName);
	// Here we check whether the function exists already or not
	// 1. find the lower bound of the map
	SurfaceFunctionMap::iterator iter_chg = functions.lower_bound(chgFuncName);

    // 2. if iter_chg == end, or if iter_chg is not a match,
    //    then this element was not in the map, so we need to insert it
	if ( iter_chg == functions.end()  ||  iter_chg->first != chgFuncName )
	{// move iter_chg to the element preceeding the insertion point
	    if ( iter_chg != functions.begin() )
	        --iter_chg;

	    // insert it
		SurfaceFunction * func = new SurfaceFunction(chgFuncName, _cavity->size());
		SurfaceFunctionPair insertion = SurfaceFunctionMap::value_type(chgFuncName, func);
		iter_chg = functions.insert(iter_chg, insertion);
    }
	
	// If it already exists there's no problem, we will pass a reference to its values to
	// _solver->compCharge(const Eigen::VectorXd &, Eigen::VectorXd &) so they will be automagically updated!
	
	_solver->compCharge(iter_pot->second->getVector(), iter_chg->second->getVector());
	}

// Revise this function. It's just a dirty hack now.
extern "C" void comp_pol_ene_pcm_(double * energy, int * separate_or_total) 
{
        if (*separate_or_total == 0) 
	{ // Using separate potentials and charges
		SurfaceFunctionMap::const_iterator iter_nuc_pot = functions.find("NucMEP");
		SurfaceFunctionMap::const_iterator iter_nuc_chg = functions.find("NucASC");
		SurfaceFunctionMap::const_iterator iter_ele_pot = functions.find("EleMEP");
		SurfaceFunctionMap::const_iterator iter_ele_chg = functions.find("EleASC");

		double UNN = (*iter_nuc_pot->second) *  (*iter_nuc_chg->second);
		double UEN = (*iter_ele_pot->second) *  (*iter_nuc_chg->second);
		double UNE = (*iter_nuc_pot->second) *  (*iter_ele_chg->second);
		double UEE = (*iter_ele_pot->second) *  (*iter_ele_chg->second);
		
		printf("U_ee = %.10E, U_en = %.10E, U_ne = %.10E, U_nn = %.10E\n", UEE, UEN, UNE, UNN);

		*energy = 0.5 * ( UNN + UEN + UNE + UEE );
        } 
	else 
	{
		SurfaceFunctionMap::const_iterator iter_pot = functions.find("TotMEP");
		SurfaceFunctionMap::const_iterator iter_chg = functions.find("TotASC");
	
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

extern "C" void get_cavity_size_(int * nts) 
{
	*nts = _cavity->size();
}

extern "C" void get_tess_centers_(double * centers) 
{
	int j = 0;
	for (int i = 0; i < _cavity->size(); i++) 
	{
		Eigen::Vector3d tess = _cavity->getElementCenter(i);
		// centers is a vector with size (3 * nr_points)
		// in which we store our tesserae centers COLUMN BY COLUMN
		// i.e. column-major order as used in Fortran.
		// The shape of the matrix in Fortran is thus centers(3, nr_points)
		centers[j] = tess(0);
		centers[j+1] = tess(1);
		centers[j+2] = tess(2);
		j += 3;
	}
	
}

extern "C" void get_tess_cent_coord_(int * its, double * center) 
{
	Eigen::Vector3d tess = _cavity->getElementCenter(*its-1);
	std::cout << tess.transpose() << std::endl;
	center[0] = tess(0);
	center[1] = tess(1);
	center[2] = tess(2);
}

extern "C" void print_pcm_()
{
	// I don't think this will work with wavelets as of now (8/7/13)
	// we should work towards this though: "Program to an interface, not an implementation."
	std::cout << "~~~~~~~~~~ PCMSolver ~~~~~~~~~~" << std::endl;
	std::cout << "========== Cavity " << std::endl;
	std::cout << *_cavity << std::endl;
	std::cout << "========== Solver " << std::endl;
	std::cout << *_solver << std::endl;
	std::cout << "============ Medium " << std::endl;
	bool fromSolvent = Input::TheInput().fromSolvent();
	if (fromSolvent)
	{
		std::cout << "Medium initialized from solvent built-in data." << std::endl;
		Solvent solvent = Input::TheInput().getSolvent();
		std::cout << solvent << std::endl;
	}
	std::cout << ".... Inside " << std::endl;
	std::cout << *(_solver->getGreenInside()) << std::endl;
	std::cout << ".... Outside " << std::endl;
	std::cout << *(_solver->getGreenOutside()) << std::endl;
	std::cout << std::endl;
}

extern "C" void print_gepol_cavity_()
{
	cout << "Cavity size" << _cavity->size() << endl;
}

extern "C" void set_surface_function_(int * nts, double * values, char * name)
{
	int nTess = _cavity->size();
	if ( nTess != *nts )
		throw std::runtime_error("You are trying to allocate a SurfaceFunction bigger than the cavity!");

	std::string functionName(name);

	// Here we check whether the function exists already or not
	// 1. find the lower bound of the map
	SurfaceFunctionMap::iterator iter = functions.lower_bound(functionName);

    // 2. if iter == end, or if iter is not a match, 
    //    then this element was not in the map, so we need to insert it
	if ( iter == functions.end()  ||  iter->first != functionName )
	{// move iter to the element preceeding the insertion point
	    if ( iter != functions.begin() )
	        --iter;

	    // insert it
		SurfaceFunction * func = new SurfaceFunction(functionName, *nts, values);
		SurfaceFunctionPair insertion = SurfaceFunctionMap::value_type(functionName, func);
		iter = functions.insert(iter, insertion);
    }
    else
    {
    	iter->second->setValues(values);
    }
}

extern "C" void get_surface_function_(int * nts, double * values, char * name) 
{
    	int nTess = _cavity->size();
	if ( nTess != *nts ) 
		throw std::runtime_error("You are trying to access a SurfaceFunction bigger than the cavity!");
	
	std::string functionName(name);
	
	SurfaceFunctionMap::const_iterator iter = functions.find(functionName);
	if ( iter == functions.end() )
		throw std::runtime_error("You are trying to access a non-existing SurfaceFunction.");

	for ( int i = 0; i < nTess; ++i )
	{
		values[i] = iter->second->getValue(i); 
	}
}

extern "C" void add_surface_function_(char * result, double * coeff, char * part) 
{
	std::string resultName(result);
	std::string partName(part);

	append_surf_func_(result);
	
	SurfaceFunctionMap::const_iterator iter_part = functions.find(partName);
	SurfaceFunctionMap::const_iterator iter_result = functions.find(resultName);

	// Using iterators and operator overloading: so neat!
	(*iter_result->second) += (*coeff) * (*iter_part->second);
}

extern "C" void print_surface_function_(char * name) 
{
	std::string functionName(name);

	SurfaceFunctionMap::const_iterator iter = functions.find(name);

	std::cout << *(iter->second) << std::endl;
}

extern "C" bool surf_func_exists_(char * name) 
{
	std::string functionName(name);

	SurfaceFunctionMap::const_iterator iter = functions.find(name);

	return iter != functions.end();
}

extern "C" void clear_surf_func_(char* name) 
{
	std::string functionName(name);

	SurfaceFunctionMap::const_iterator iter = functions.find(name);

	iter->second->clear();
}

extern "C" void append_surf_func_(char* name) 
{
	int nTess = _cavity->size();
	std::string functionName(name);

	// Here we check whether the function exists already or not
	// 1. find the lower bound of the map
	SurfaceFunctionMap::iterator iter = functions.lower_bound(functionName);

    // 2. if iter == end, or if iter is not a match, 
    //    then this element was not in the map, so we need to insert it
	if ( iter == functions.end()  ||  iter->first != functionName )
	{// move iter to the element preceeding the insertion point
	    if ( iter != functions.begin() )
	        --iter;

	    // insert it
		SurfaceFunction * func = new SurfaceFunction(functionName, nTess);
		SurfaceFunctionPair insertion = SurfaceFunctionMap::value_type(functionName, func);
		iter = functions.insert(iter, insertion);
    }
    else
    {// What happens if it is already in the map? The values need to be updated.
     // Nothing, I assume that if one calls append_surf_func_ will then also call
     // set_surface_function_ somewhere else, hence the update will be done there.
 	} 
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
	std::vector<Sphere> spheres;

	if (_mode == "Implicit") 
	{
		initSpheresImplicit(charges, centers, spheres);
		parsedInput.setSpheres(spheres);
	} 
	else if (_mode == "Atoms") 
	{
		initSpheresAtoms(charges, centers, spheres);
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
	// TODO: since WaveletCavity extends cavity in a significant way, use of the Factory Method design pattern does not work for wavelet cavities. (8/7/13)
	std::string modelType = Input::TheInput().getSolverType();
        if (modelType == "Wavelet" || modelType == "Linear") 
	{// Both PWC and PWL require a WaveletCavity
		_waveletCavity = initWaveletCavity();
	}
        else
	{// This means in practice that the CavityFactory is now working only for GePol.
		_cavity = CavityFactory::TheCavityFactory().createCavity(cavityType, spheres, area, probeRadius, addSpheres, patchLevel, coarsity);
	}	
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
	// This thing is rather ugly I admit, but will be changed (as soon as wavelet PCM is working with DALTON)
	// it is needed because: 1. comment above on cavities; 2. wavelet cavity and solver depends on each other
	// (...not our fault, but should remedy somehow)
        if (modelType == "Wavelet") 
	{
		_PWCSolver = new PWCSolver(gfInside, gfOutside);
		_PWCSolver->buildSystemMatrix(*_waveletCavity);
		_waveletCavity->uploadPoints(_PWCSolver->getQuadratureLevel(), _PWCSolver->getT_(), false); // WTF is happening here???
		_cavity = _waveletCavity;
		_solver = _PWCSolver;
	} 
	else if (modelType == "Linear") 
	{
		_PWLSolver = new PWLSolver(gfInside, gfOutside);
		_PWLSolver->buildSystemMatrix(*_waveletCavity);
		_waveletCavity->uploadPoints(_PWLSolver->getQuadratureLevel(),_PWLSolver->getT_(), true); // WTF is happening here???
		_cavity = _waveletCavity;
		_solver = _PWLSolver;
	}
        else
	{// This means that the factory is properly working only for IEFSolver and CPCMSolver
		_solver = SolverFactory::TheSolverFactory().createSolver(modelType, gfInside, gfOutside, correction, eqType);
		_solver->buildSystemMatrix(*_cavity);
	}	
}

void initAtoms(Eigen::VectorXd & charges_, Eigen::Matrix3Xd & sphereCenter_) 
{
	int nuclei;
	collect_nctot_(&nuclei);
	sphereCenter_.resize(Eigen::NoChange, nuclei);
	charges_.resize(nuclei);
	double * chg = charges_.data();
	double * centers = sphereCenter_.data();
	collect_atoms_(chg, centers);
} 

void initSpheresAtoms(const Eigen::VectorXd & charges_, const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_) 
{
	vector<int> atomsInput = Input::TheInput().getAtoms();
	vector<double> radiiInput = Input::TheInput().getRadii();
	
	for (int i = 0; i < atomsInput.size(); ++i) 
	{
		int index = atomsInput[i] - 1; // -1 to go from human readable to machine readable
		Eigen::Vector3d center = sphereCenter_.col(index);
		Sphere sph(center, radiiInput[i]);
		spheres_.push_back(sph);
	}
}

void initSpheresImplicit(const Eigen::VectorXd & charges_, const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_) 
{
	bool scaling = Input::TheInput().getScaling();
	std::string set = Input::TheInput().getRadiiSet();
	
	std::vector<Atom> radiiSet;
	if ( set == "UFF" )
	{
		radiiSet = Atom::initUFF();
	}
	else
	{
		radiiSet = Atom::initBondi();
	}

	for (int i = 0; i < charges_.size(); ++i) 
	{
		int index = charges_(i) - 1;
		double radius = radiiSet[index].getAtomRadius();
                if (scaling) 
		{
			radius *= radiiSet[index].getAtomRadiusScaling();
                }
		Eigen::Vector3d center = sphereCenter_.col(i);
		Sphere sph(center, radius);
		spheres_.push_back(sph);
	}
}

WaveletCavity * initWaveletCavity()
{
	int patchLevel = Input::TheInput().getPatchLevel();
	std::vector<Sphere> spheres = Input::TheInput().getSpheres();
	double coarsity = Input::TheInput().getCoarsity();
	double probeRadius = Input::TheInput().getProbeRadius();
    	
	WaveletCavity * cav = new WaveletCavity(spheres, probeRadius, patchLevel, coarsity);
	cav->readCavity("molec_dyadic.dat");
	
	return cav;
}