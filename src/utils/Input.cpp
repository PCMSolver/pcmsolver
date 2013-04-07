#include "Input.h"

Input::Input(const char *parsedInputFile)
{
	// Create a Getkw object from input file.
        Getkw input = Getkw(parsedInputFile, false, true);
 	
	const Section & cavity = input.getSect("Cavity<Cavity>");

        type = cavity.getStr("Type");
	if (type == "GePol") {
		area = cavity.getDbl("Area");
	} else {
		patchLevel = cavity.getDbl("PatchLevel"); 
		coarsity = cavity.getDbl("Coarsity");
	}
	scaling = cavity.getBool("Scaling");
	addSpheres = cavity.getBool("AddSpheres");
	mode = cavity.getStr("Mode");
 	atoms = cavity.getIntVec("Atoms");
	radii = cavity.getDblVec("Radii");
	// spheres must be treated with some special care
	//std::vector<Sphere> spheres;
        
	const Section & medium = input.getSect("Medium<Medium>");

	solvent = medium.getStr("Solvent");
	solverType = medium.getStr("SolverType");
	equationType = medium.getStr("EquationType");
	probeRadius = medium.getDbl("ProbeRadius");

		/*std::string greenType;
		std::string derivativeType;
		double epsilon;
		double epsilonReal;
		double epsilonImaginary;
		std::vector<double> spherePosition;
		double sphereRadius;*/

}
