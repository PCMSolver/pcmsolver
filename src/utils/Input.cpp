#include "Input.hpp"

#include <map>
#include <vector>
#include <stdexcept>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif
#include "Getkw.h"

#include "Solvent.hpp"
#include "Sphere.hpp"

Input::Input()
{
	const char * parsedInputFile = "@pcmsolver.inp";
	// Create a Getkw object from input file.
        Getkw input = Getkw(parsedInputFile, false, true);
	const Section & cavity = input.getSect("Cavity");

	cavFilename = cavity.getStr("Restart");

        type = cavity.getStr("Type");
	if (type == "GePol") 
	{ // GePol cavity branch
		area = cavity.getDbl("Area");
		patchLevel = 0.0;
		coarsity = 0.0;
	}
	else if (type == "Wavelet") 
	{ // Wavelet cavity branch
#if defined (WAVELET_DEVELOPMENT)
		area = 0.0;
		patchLevel = cavity.getInt("PatchLevel"); 
		coarsity = cavity.getDbl("Coarsity");
#else
		throw std::runtime_error("Wavelet cavity generator is not included in this release.");
#endif
	}
	else
	{ // TsLess cavity branch
#if defined (TSLESS_DEVELOPMENT)
		// Get the necessary input parameters to be passed to the TsLessCavity CTOR.
#else
		throw std::runtime_error("TsLess cavity generator is not included in this release.");
#endif
	}
	scaling = cavity.getBool("Scaling");
	radiiSet = cavity.getStr("RadiiSet");
	addSpheres = cavity.getBool("AddSpheres");
	mode = cavity.getStr("Mode");
	if (mode == "Explicit") 
	{
		std::vector<double> spheresInput = cavity.getDblVec("Spheres");
		int j = 0;
		int upperBound = (int)spheresInput.size() / 4;
		for (int i = 0; i < upperBound; ++i) {
			Eigen::Vector3d center;
			center << spheresInput[j], spheresInput[j+1], spheresInput[j+2];
			Sphere sph(center, spheresInput[j+3]);
			spheres.push_back(sph);
			j += 4;
		}
	} 
	else if (mode == "Atoms")
	{
 		atoms = cavity.getIntVec("Atoms");
		radii = cavity.getDblVec("Radii");
	}
       
        // Get the contents of the Medium section	
	const Section & medium = input.getSect("Medium");
	// Get the name of the solvent
	std::string _name = medium.getStr("Solvent");
	// A useful map to convert the Der string to an integer
	// which will be passed to the Green's function CTOR.
	std::map<std::string, int> mapStringToIntDerType;
	mapStringToIntDerType.insert(std::map<std::string, int>::value_type("Numerical", 0));
	mapStringToIntDerType.insert(std::map<std::string, int>::value_type("Analytic", 1));
	mapStringToIntDerType.insert(std::map<std::string, int>::value_type("Derivative", 2));
	mapStringToIntDerType.insert(std::map<std::string, int>::value_type("Gradient", 3));
	mapStringToIntDerType.insert(std::map<std::string, int>::value_type("Hessian", 4));

	if (_name == "Explicit") 
	{
		hasSolvent = false;
		// Get the probe radius
		probeRadius = medium.getDbl("ProbeRadius");
		// Get the contents of the Green<inside> section...
		const Section & _inside = medium.getSect("Green<inside>");
		// ...and initialize the data members
		greenInsideType = _inside.getStr("Type");
		derivativeInsideType = mapStringToIntDerType.find(_inside.getStr("Der"))->second;
		epsilonInside = _inside.getDbl("Eps");
		// Get the contents of the Green<outside> section...
		const Section & _outside = medium.getSect("Green<outside>");
		// ...and initialize the data members
		greenOutsideType = _outside.getStr("Type");
		derivativeOutsideType = mapStringToIntDerType.find(_outside.getStr("Der"))->second;
		epsilonOutside = _outside.getDbl("Eps");
		// This will be needed for the metal sphere
		epsilonReal = medium.getDbl("EpsRe");
		epsilonImaginary = medium.getDbl("EpsImg");
		spherePosition = medium.getDblVec("SpherePosition");
		sphereRadius = medium.getDbl("SphereRadius");
	} 
	else // This part must be reviewed!! Some data members are not initialized... 
	{       // Just initialize the solvent object in this class
		hasSolvent = true;
		std::map<std::string, Solvent> solvents = Solvent::initSolventMap();
		solvent = solvents[_name];
		probeRadius = solvent.getRadius();
		// Specification of the solvent by name means isotropic PCM
		// We have to initialize the Green's functions data here, Solvent class
		// is an helper class and should not be used in the core classes.
		greenInsideType = "Vacuum";
		derivativeInsideType = mapStringToIntDerType.find("Derivative")->second;
		epsilonInside = 1.0;
	        greenOutsideType = "UniformDielectric";
		derivativeOutsideType = mapStringToIntDerType.find("Derivative")->second;	
		epsilonOutside = solvent.getEpsStatic();
	}
	// A useful map to convert the EquationType string to an integer
	// which will be passed to the Solver CTOR.
	std::map<std::string, int> mapStringToIntEquationType;
	mapStringToIntEquationType.insert(std::map<std::string, int>::value_type("FirstKind", 0));
	mapStringToIntEquationType.insert(std::map<std::string, int>::value_type("SecondKind", 1));
	mapStringToIntEquationType.insert(std::map<std::string, int>::value_type("Full", 2));

	solverType = medium.getStr("SolverType");
	equationType = mapStringToIntEquationType.find(medium.getStr("EquationType"))->second;
	correction = medium.getDbl("Correction");

#if not defined (WAVELET_DEVELOPMENT)
	if (solverType == "Wavelet" || solverType == "Linear")
	{
		throw std::runtime_error("Wavelet cavity generator and solver are not included in this release.");
	}
#endif

	// Now we have all input parameters, do some sanity checks
        if ( type == "GePol" || type == "TsLess" )
	{
		if (solverType == "Wavelet" || solverType == "Linear") // User asked for GePol or TsLess cavity with wavelet solvers
		{
	    		throw std::runtime_error("GePol cavity can be used only with traditional solvers.");
		}
	}
	else 
	{
		if (solverType == "IEFPCM" || solverType == "CPCM") // User asked for wavelet cavity with traditional solvers
		{
	    		throw std::runtime_error("Wavelet cavity can be used only with wavelet solvers.");
		}
	}
}

