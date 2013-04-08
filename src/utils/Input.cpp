#include "Input.h"

Input::Input(const char *parsedInputFile)
{
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
	// spheres must be treated with some special care. Maybe move some interface level functions here?
	//std::vector<Sphere> spheres;
        
	const Section & medium = input.getSect("Medium");
	// Create a Solvent object here?
	solvent = medium.getStr("Solvent");
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
