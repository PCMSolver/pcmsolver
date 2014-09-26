/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *     
 *     This file is part of PCMSolver.
 *     
 *     PCMSolver is free software: you can redistribute it and/or modify       
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *     
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *     
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *     
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#include "Input.hpp"

#include <map>
#include <vector>
#include <stdexcept>

#include "Config.hpp"

#include <Eigen/Dense>
#include "Getkw.h"

#include <boost/algorithm/string.hpp>

#include "InputManager.hpp"
#include "PhysicalConstants.hpp"
#include "Solvent.hpp"
#include "Sphere.hpp"

using boost::algorithm::to_upper_copy;

Input::Input(const char * pythonParsed)
{
    reader(pythonParsed);
    semanticParser();
}

Input::Input(const cavityInput & cav, const solverInput & solv, const greenInput & green)
{
    reader(cav, solv, green);
    semanticParser();
}

void Input::reader(const char * pythonParsed)
{
    // Create a Getkw object from input file.
    Getkw input = Getkw(pythonParsed, false, true);

    CODATAyear_ = input.getInt("CODATA");

    const Section & cavity = input.getSect("Cavity");

    type_ = cavity.getStr("Type");
    area_ = cavity.getDbl("Area");
    patchLevel_ = cavity.getInt("PatchLevel");
    coarsity_ = cavity.getDbl("Coarsity");
    minDistance_ = cavity.getDbl("MinDistance");
    derOrder_ = cavity.getInt("DerOrder");
    cavFilename_ = cavity.getStr("Restart");

    scaling_ = cavity.getBool("Scaling");
    radiiSet_ = cavity.getStr("RadiiSet");
    minimalRadius_ = cavity.getDbl("MinRadius");
    mode_ = cavity.getStr("Mode");
    if (mode_ == "Explicit") {
        std::vector<double> spheresInput = cavity.getDblVec("Spheres");
        int j = 0;
        int upperBound = int(spheresInput.size() / 4);
        for (int i = 0; i < upperBound; ++i) {
            Eigen::Vector3d center;
            center << spheresInput[j], spheresInput[j+1], spheresInput[j+2];
            Sphere sph(center, spheresInput[j+3]);
            spheres_.push_back(sph);
            j += 4;
        }
    } else if (mode_ == "Atoms") {
        atoms_ = cavity.getIntVec("Atoms");
        radii_ = cavity.getDblVec("Radii");
    }

    // Get the contents of the Medium section
    const Section & medium = input.getSect("Medium");
    // Get the name of the solvent
    std::string name = medium.getStr("Solvent");
    if (name == "Explicit") {
        hasSolvent_ = false;
        // Get the probe radius
        probeRadius_ = medium.getDbl("ProbeRadius");
        // Get the contents of the Green<inside> section...
        const Section & inside = medium.getSect("Green<inside>");
        // ...and initialize the data members
        greenInsideType_ = inside.getStr("Type");
	derivativeInsideType_ = derivativeTraits(inside.getStr("Der"));
	epsilonInside_ = inside.getDbl("Eps");
        // Get the contents of the Green<outside> section...
        const Section & outside = medium.getSect("Green<outside>");
        // ...and initialize the data members
        greenOutsideType_ = outside.getStr("Type");
        derivativeOutsideType_ = derivativeTraits(inside.getStr("Der"));
        epsilonOutside_ = outside.getDbl("Eps");
        // This will be needed for the metal sphere only
	if (greenOutsideType_ == "MetalSphere") {
	        epsilonReal_ = outside.getDbl("EpsRe");
        	epsilonImaginary_ = outside.getDbl("EpsImg");
	        spherePosition_ = outside.getDblVec("SpherePosition");
        	sphereRadius_ = outside.getDbl("SphereRadius");
	}
    } else { // This part must be reviewed!! Some data members are not initialized...
        // Just initialize the solvent object in this class
        hasSolvent_ = true;
        std::map<std::string, Solvent> solvents = Solvent::initSolventMap();
        solvent_ = solvents[name];
        probeRadius_ = solvent_.probeRadius() * angstromToBohr(CODATAyear_);
        // Specification of the solvent by name means isotropic PCM
        // We have to initialize the Green's functions data here, Solvent class
        // is an helper class and should not be used in the core classes.
        greenInsideType_ = "Vacuum";
        derivativeInsideType_ = derivativeTraits("Derivative");
        epsilonInside_ = 1.0;
        greenOutsideType_ = "UniformDielectric";
        derivativeOutsideType_ = derivativeTraits("Derivative");
        epsilonOutside_ = solvent_.epsStatic();
    }
    
    solverType_ = medium.getStr("SolverType");
    equationType_ = integralEquation(medium.getStr("EquationType"));
    correction_ = medium.getDbl("Correction");
    hermitivitize_ = medium.getBool("MatrixSymm");
    
    providedBy_ = std::string("API-side");
}

void Input::reader(const cavityInput & cav, const solverInput & solv, const greenInput & green)
{
    CODATAyear_ = 2010;

    type_ = to_upper_copy(std::string(cav.cavity_type));
    area_ = cav.area;
    patchLevel_ = cav.patch_level;
    coarsity_ = cav.coarsity;
    minDistance_ = cav.min_distance;
    derOrder_ = cav.der_order;
    cavFilename_ = std::string(cav.restart_name); // No case conversion here!

    scaling_ = cav.scaling;
    radiiSet_ = cav.radii_set;
    minimalRadius_ = cav.min_radius;
    mode_ = std::string("IMPLICIT"); 

    std::string name = to_upper_copy(std::string(solv.solvent));
    if (name.empty()) {
        hasSolvent_ = false;
        // Get the probe radius
        probeRadius_ = solv.probe_radius;
        // Get the contents of the Green<inside> section...
        // ...and initialize the data members
        greenInsideType_ = to_upper_copy(std::string(green.inside_type));
	derivativeInsideType_ = derivativeTraits("DERIVATIVE");
	epsilonInside_ = 1.0;
        // Get the contents of the Green<outside> section...
        // ...and initialize the data members
        greenOutsideType_ = to_upper_copy(std::string(green.outside_type));
        derivativeOutsideType_ = derivativeTraits("DERIVATIVE");
        epsilonOutside_ = green.outside_epsilon;
    } else { // This part must be reviewed!! Some data members are not initialized...
        // Just initialize the solvent object in this class
        hasSolvent_ = true;
        std::map<std::string, Solvent> solvents = Solvent::initSolventMap();
        solvent_ = solvents[name];
        probeRadius_ = solvent_.probeRadius() * angstromToBohr(CODATAyear_);
        // Specification of the solvent by name means isotropic PCM
        // We have to initialize the Green's functions data here, Solvent class
        // is an helper class and should not be used in the core classes.
        greenInsideType_ = std::string("VACUUM");
        derivativeInsideType_ = derivativeTraits("DERIVATIVE");
        epsilonInside_ = 1.0;
        greenOutsideType_ = std::string("UNIFORMDIELECTRIC");
        derivativeOutsideType_ = derivativeTraits("DERIVATIVE");
        epsilonOutside_ = solvent_.epsStatic();
    }
    
    
    solverType_ = to_upper_copy(std::string(solv.solver_type));
    std::string inteq = to_upper_copy(std::string(solv.equation_type));
    equationType_ = integralEquation(inteq);
    correction_ = solv.correction;
    hermitivitize_ = solv.matrix_symmetrize;
    
    providedBy_ = std::string("host-side");
}

void Input::semanticParser()
{
#if not defined (WAVELET_DEVELOPMENT) || not defined (TSLESS_DEVELOPMENT)
    if (type_ == "GEPOL" || type_ == "TSLESS") {
        if (solverType_ == "WAVELET" || solverType_ == "LINEAR") { // User asked for GePol or TsLess cavity with wavelet solvers
            throw std::runtime_error("GePol cavity can be used only with traditional solvers.");
        }
    } else if (type_ == "WAVELET" ) { // Hoping that the user knows what's going on if he asks for a restart...
        if (solverType_ == "IEFPCM" || solverType_ == "CPCM") { // User asked for wavelet cavity with traditional solvers
            throw std::runtime_error("Wavelet cavity can be used only with wavelet solvers.");
        }
    }
#endif
}

int derivativeTraits(const std::string & name)
{
    static std::map<std::string, int> mapStringToIntDerType;
    mapStringToIntDerType.insert(std::map<std::string, int>::value_type("NUMERICAL", 0));
    mapStringToIntDerType.insert(std::map<std::string, int>::value_type("DERIVATIVE", 1));
    mapStringToIntDerType.insert(std::map<std::string, int>::value_type("GRADIENT", 2));
    mapStringToIntDerType.insert(std::map<std::string, int>::value_type("HESSIAN", 3));

    return mapStringToIntDerType.find(name)->second;
}

int integralEquation(const std::string & name)
{
    static std::map<std::string, int> mapStringToIntEquationType;
    mapStringToIntEquationType.insert(std::map<std::string, int>::value_type("FIRSTKIND",
                                      0));
    mapStringToIntEquationType.insert(
        std::map<std::string, int>::value_type("SECONDKIND", 1));
    mapStringToIntEquationType.insert(std::map<std::string, int>::value_type("FULL", 2));
    
    return mapStringToIntEquationType.find(name)->second;
}
