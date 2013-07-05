#include "Input.hpp"

#include <map>

#include "Getkw.h"

#include "Solvent.hpp"
#include "Sphere.hpp"

Input::Input()
{
	const char * parsedInputFile = "@pcmsolver.inp";
	// Create a Getkw object from input file.
        Getkw input = Getkw(parsedInputFile, false, true);
 	
	const Section & cavity = input.getSect("Cavity");

        type = cavity.getStr("Type");
	if (type == "GePol") { 
		area = cavity.getDbl("Area");
		patchLevel = 0.0;
		coarsity = 0.0;
	} else {
		area = 0.0;
		patchLevel = cavity.getDbl("PatchLevel"); 
		coarsity = cavity.getDbl("Coarsity");
	}
	scaling = cavity.getBool("Scaling");
	addSpheres = cavity.getBool("AddSpheres");
	mode = cavity.getStr("Mode");
	if (mode == "Explicit") 
	{
		std::vector<double> spheresInput = cavity.getDblVec("Spheres");
		int j = 0;
		for (int i = 0; i < spheresInput.size(); ++i) {
			Eigen::Vector3d center;
			center << spheresInput[j], spheresInput[j+1], spheresInput[j+2];
			Sphere sph(center, spheresInput[j+3]);
			spheres.push_back(sph);
			j += 4;
		}
	} else if (mode == "Atoms")
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
}

