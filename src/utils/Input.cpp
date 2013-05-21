#include "Input.h"

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
 	atoms = cavity.getIntVec("Atoms");
	radii = cavity.getDblVec("Radii");
	if (mode == "Explicit") {
		std::vector<double> spheresInput = cavity.getDblVec("Spheres");
		int j = 0;
		for (int i = 0; i < spheresInput.size(); ++i) {
			Eigen::Vector3d center;
			center << spheresInput[j], spheresInput[j+1], spheresInput[j+2];
			Sphere sph(center, spheresInput[j+3]);
			spheres.push_back(sph);
			j += 4;
		}
	}
        
	const Section & medium = input.getSect("Medium");
	// Create a Solvent object here?
	std::string _name = medium.getStr("Solvent");
	SolventMap solvents = Solvent::initSolventMap();
	solvent = solvents[_name];

	solverType = medium.getStr("SolverType");
	equationType = medium.getStr("EquationType");
	correction = medium.getDbl("Correction");
	probeRadius = medium.getDbl("ProbeRadius");
	greenType = medium.getStr("Type");
	derivativeType = medium.getStr("Der");
	epsilon = medium.getDbl("Eps");
	epsilonReal = medium.getDbl("EpsRe");
	epsilonImaginary = medium.getDbl("EpsImg");
	spherePosition = medium.getDblVec("SpherePosition");
	sphereRadius = medium.getDbl("SphereRadius");
}

