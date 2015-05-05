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

#include <iostream>
#include <map>
#include <vector>
#include <stdexcept>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>
#include "Getkw.h"

#include <boost/algorithm/string.hpp>

#include "CavityData.hpp"
#include "GreenData.hpp"
#include "SolverData.hpp"
#include "InputManager.hpp"
#include "PhysicalConstants.hpp"
#include "Solvent.hpp"
#include "Sphere.hpp"

using boost::algorithm::to_upper_copy;

Input::Input(const std::string & filename)
{
    reader(filename.c_str());
    semanticCheck();
}

Input::Input(const cavityInput & cav, const solverInput & solv,
             const greenInput & green)
{
    reader(cav, solv, green);
    semanticCheck();
}

void Input::reader(const char * pythonParsed)
{
    // Create a Getkw object from input file.
    Getkw input = Getkw(pythonParsed, false, true);

    units_      = input.getStr("UNITS");
    CODATAyear_ = input.getInt("CODATA");

    const Section & cavity = input.getSect("CAVITY");

    type_ = cavity.getStr("TYPE");
    area_ = cavity.getDbl("AREA");
    patchLevel_ = cavity.getInt("PATCHLEVEL");
    coarsity_ = cavity.getDbl("COARSITY");
    minDistance_ = cavity.getDbl("MINDISTANCE");
    derOrder_ = cavity.getInt("DERORDER");
    if (type_ == "RESTART") {
        cavFilename_ = cavity.getStr("NPZFILE");
    }

    scaling_ = cavity.getBool("SCALING");
    radiiSet_ = to_upper_copy(cavity.getStr("RADIISET"));
    minimalRadius_ = cavity.getDbl("MINRADIUS");
    mode_ = to_upper_copy(cavity.getStr("MODE"));
    if (mode_ == "EXPLICIT") {
        std::vector<double> spheresInput = cavity.getDblVec("SPHERES");
        int j = 0;
        int nAtoms = int(spheresInput.size() / 4);
        for (int i = 0; i < nAtoms; ++i) {
            Eigen::Vector3d center;
            center << spheresInput[j], spheresInput[j+1], spheresInput[j+2];
            Sphere sph(center, spheresInput[j+3]);
            spheres_.push_back(sph);
            j += 4;
        }
        molecule_ = Molecule(spheres_);
    } else if (mode_ == "ATOMS") {
        atoms_ = cavity.getIntVec("ATOMS");
        radii_ = cavity.getDblVec("RADII");
    }

    // Get the contents of the Medium section
    const Section & medium = input.getSect("MEDIUM");
    // Get the name of the solvent
    std::string name = medium.getStr("SOLVENT");
    if (name == "EXPLICIT") {
        hasSolvent_ = false;
        // Get the probe radius
        probeRadius_ = medium.getDbl("PROBERADIUS");
        // Get the contents of the Green<inside> section...
        const Section & inside = medium.getSect("GREEN<INSIDE>");
        // ...and initialize the data members
        greenInsideType_ = inside.getStr("TYPE");
        derivativeInsideType_ = derivativeTraits(inside.getStr("DER"));
        epsilonInside_ = inside.getDbl("EPS");
        // Get the contents of the Green<outside> section...
        const Section & outside = medium.getSect("GREEN<OUTSIDE>");
        // ...and initialize the data members
        greenOutsideType_ = outside.getStr("TYPE");
        derivativeOutsideType_ = derivativeTraits(outside.getStr("DER"));
        epsilonStaticOutside_ = outside.getDbl("EPS");
        epsilonDynamicOutside_ = outside.getDbl("EPSDYN");
        // This will be needed for the metal sphere only
        if (greenOutsideType_ == "METALSPHERE") {
            epsilonReal_ = outside.getDbl("EPSRE");
            epsilonImaginary_ = outside.getDbl("EPSIMG");
            spherePosition_ = outside.getDblVec("SPHEREPOSITION");
            sphereRadius_ = outside.getDbl("SPHERERADIUS");
        }
    } else { // This part must be reviewed!! Some data members are not initialized...
        // Just initialize the solvent object in this class
        hasSolvent_ = true;
        std::map<std::string, Solvent> solvents = Solvent::initSolventMap();
        if (solvents.find(name) == solvents.end()) {
            throw std::runtime_error("Solvent " + name + " NOT found!");
        } else {
            solvent_ = solvents[name];
        }
        probeRadius_ = solvent_.probeRadius() * angstromToBohr(CODATAyear_);
        // Specification of the solvent by name means isotropic PCM
        // We have to initialize the Green's functions data here, Solvent class
        // is an helper class and should not be used in the core classes.
        greenInsideType_ = "VACUUM";
        derivativeInsideType_ = derivativeTraits("DERIVATIVE");
        epsilonInside_ = 1.0;
        greenOutsideType_ = "UNIFORMDIELECTRIC";
        derivativeOutsideType_ = derivativeTraits("DERIVATIVE");
        epsilonStaticOutside_ = solvent_.epsStatic();
        epsilonDynamicOutside_ = solvent_.epsDynamic();
    }
    integratorType_ = "COLLOCATION"; // Currently hardcoded!!!

    solverType_ = medium.getStr("SOLVERTYPE");
    equationType_ = integralEquation(medium.getStr("EQUATIONTYPE"));
    correction_ = medium.getDbl("CORRECTION");
    hermitivitize_ = medium.getBool("MATRIXSYMM");

    providedBy_ = std::string("API-side");
}

void Input::reader(const cavityInput & cav, const solverInput & solv,
                   const greenInput & green)
{
    CODATAyear_ = 2010;

    type_ = to_upper_copy(std::string(cav.cavity_type));
    area_ = cav.area * angstrom2ToBohr2(CODATAyear_);
    patchLevel_ = cav.patch_level;
    coarsity_ = cav.coarsity * angstromToBohr(CODATAyear_);
    minDistance_ = cav.min_distance * angstromToBohr(CODATAyear_);
    derOrder_ = cav.der_order;
    if (type_ == "RESTART") {
        cavFilename_ = std::string(cav.restart_name); // No case conversion here!
    }

    scaling_ = cav.scaling;
    radiiSet_ = cav.radii_set;
    minimalRadius_ = cav.min_radius * angstromToBohr(CODATAyear_);
    mode_ = std::string("IMPLICIT");

    std::string name = to_upper_copy(std::string(solv.solvent));
    if (name.empty()) {
        hasSolvent_ = false;
        // Get the probe radius
        probeRadius_ = solv.probe_radius * angstromToBohr(CODATAyear_);
        // Get the contents of the Green<inside> section...
        // ...and initialize the data members
        greenInsideType_ = to_upper_copy(std::string(green.inside_type));
        derivativeInsideType_ = derivativeTraits("DERIVATIVE");
        epsilonInside_ = 1.0;
        // Get the contents of the Green<outside> section...
        // ...and initialize the data members
        greenOutsideType_ = to_upper_copy(std::string(green.outside_type));
        derivativeOutsideType_ = derivativeTraits("DERIVATIVE");
        epsilonStaticOutside_ = green.outside_epsilon;
        epsilonDynamicOutside_ = green.outside_epsilon;
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
        epsilonStaticOutside_ = solvent_.epsStatic();
        epsilonDynamicOutside_ = solvent_.epsDynamic();
    }
    integratorType_ = "COLLOCATION"; // Currently hardcoded!!!

    solverType_ = to_upper_copy(std::string(solv.solver_type));
    std::string inteq = to_upper_copy(std::string(solv.equation_type));
    equationType_ = integralEquation(inteq);
    correction_ = solv.correction;
    hermitivitize_ = true;

    providedBy_ = std::string("host-side");
}

void Input::semanticCheck()
{
}

cavityData Input::cavityParams()
{
    if (cavData_.empty) {
        cavData_ = cavityData(molecule_, area_, probeRadius_, minDistance_, derOrder_,
                              minimalRadius_,
                              patchLevel_, coarsity_, cavFilename_);
    }
    return cavData_;
}

greenData Input::insideGreenParams()
{
    if (insideGreenData_.empty) {
        insideGreenData_ = greenData(derivativeInsideType_, epsilonInside_, integratorType_);
    }
    return insideGreenData_;
}

greenData Input::outsideStaticGreenParams()
{
    if (outsideStaticGreenData_.empty) {
        outsideStaticGreenData_ = greenData(derivativeOutsideType_, epsilonStaticOutside_,
                                            integratorType_);
    }
    return outsideStaticGreenData_;
}

greenData Input::outsideDynamicGreenParams()
{
    if (outsideDynamicGreenData_.empty) {
        outsideDynamicGreenData_ = greenData(derivativeOutsideType_, epsilonDynamicOutside_,
                                             integratorType_);
    }
    return outsideDynamicGreenData_;
}

int derivativeTraits(const std::string & name)
{
    static std::map<std::string, int> mapStringToIntDerType;
    mapStringToIntDerType.insert(std::map<std::string, int>::value_type("NUMERICAL", 0));
    mapStringToIntDerType.insert(std::map<std::string, int>::value_type("DERIVATIVE",
                                 1));
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
