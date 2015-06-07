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
    dyadicFilename_ = cavity.getStr("DYADICFILE");

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
    if (name == "EXPLICIT" || name == "E") {
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
        epsilonStatic1_ = outside.getDbl("EPS1");
        epsilonDynamic1_ = outside.getDbl("EPSDYN1");
        epsilonStatic2_ = outside.getDbl("EPS2");
        epsilonDynamic2_ = outside.getDbl("EPSDYN2");
        center_ = outside.getDbl("CENTER");
        width_ = outside.getDbl("WIDTH");
        origin_ = outside.getDblVec("INTERFACEORIGIN");
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
    integratorType_ = integratorPolicy("COLLOCATION"); // Currently hardcoded!!!

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
    integratorType_ = integratorPolicy("COLLOCATION"); // Currently hardcoded!!!

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
        cavData_ = cavityData(molecule_, area_, probeRadius_, minDistance_, derOrder_, minimalRadius_,
                        patchLevel_, coarsity_, cavFilename_, dyadicFilename_);
    }
    return cavData_;
}

greenData Input::insideGreenParams()
{
    if (insideGreenData_.empty) {
        int profile = profilePolicy("UNIFORM");
        insideGreenData_ = greenData(derivativeInsideType_, integratorType_, profile, epsilonInside_);
    }
    return insideGreenData_;
}

greenData Input::outsideStaticGreenParams()
{
    if (outsideStaticGreenData_.empty) {
        int profile = profilePolicy("UNIFORM");
        outsideStaticGreenData_ = greenData(derivativeOutsideType_,
                                            integratorType_, profile,  epsilonStaticOutside_);
        if (not hasSolvent_) {
           outsideStaticGreenData_.howProfile  = profilePolicy("ONELAYERTANH");
           outsideStaticGreenData_.epsilon1 = epsilonStatic1_;
           outsideStaticGreenData_.epsilon2 = epsilonStatic2_;
           outsideStaticGreenData_.center   = center_;
           outsideStaticGreenData_.width    = width_;
           outsideStaticGreenData_.origin   << origin_[0], origin_[1], origin_[2];
        }
    }
    return outsideStaticGreenData_;
}

greenData Input::outsideDynamicGreenParams()
{
    if (outsideDynamicGreenData_.empty) {
        int profile = profilePolicy("UNIFORM");
        outsideDynamicGreenData_ = greenData(derivativeOutsideType_,
                                             integratorType_, profile, epsilonDynamicOutside_);
        if (not hasSolvent_) {
           outsideDynamicGreenData_.howProfile  = profilePolicy("ONELAYERTANH");
           outsideDynamicGreenData_.epsilon1 = epsilonDynamic1_;
           outsideDynamicGreenData_.epsilon2 = epsilonDynamic2_;
           outsideDynamicGreenData_.center   = center_;
           outsideDynamicGreenData_.width    = width_;
           outsideDynamicGreenData_.origin   << origin_[0], origin_[1], origin_[2];
        }
    }
    return outsideDynamicGreenData_;
}

solverData Input::solverParams()
{
    if (solverData_.empty) {
        solverData_ = solverData(correction_, equationType_, hermitivitize_);
    }
    return solverData_;
}

int derivativeTraits(const std::string & name)
{
    static std::map<std::string, int> mapStringToInt;
    mapStringToInt.insert(std::map<std::string, int>::value_type("NUMERICAL", 0));
    mapStringToInt.insert(std::map<std::string, int>::value_type("DERIVATIVE", 1));
    mapStringToInt.insert(std::map<std::string, int>::value_type("GRADIENT", 2));
    mapStringToInt.insert(std::map<std::string, int>::value_type("HESSIAN", 3));

    return mapStringToInt.find(name)->second;
}

int integratorPolicy(const std::string & name)
{
    static std::map<std::string, int> mapStringToInt;
    mapStringToInt.insert(std::map<std::string, int>::value_type("COLLOCATION", 0));
    mapStringToInt.insert(std::map<std::string, int>::value_type("NUMERICAL", 1));

    return mapStringToInt.find(name)->second;
}

int profilePolicy(const std::string & name)
{
    static std::map<std::string, int> mapStringToInt;
    mapStringToInt.insert(std::map<std::string, int>::value_type("UNIFORM", 0));
    mapStringToInt.insert(std::map<std::string, int>::value_type("YUKAWA", 1));
    mapStringToInt.insert(std::map<std::string, int>::value_type("ANISOTROPIC", 2));
    mapStringToInt.insert(std::map<std::string, int>::value_type("ONELAYERTANH", 3));

    return mapStringToInt.find(name)->second;
}

int integralEquation(const std::string & name)
{
    static std::map<std::string, int> mapStringToInt;
    mapStringToInt.insert(std::map<std::string, int>::value_type("FIRSTKIND",0));
    mapStringToInt.insert(std::map<std::string, int>::value_type("SECONDKIND", 1));
    mapStringToInt.insert(std::map<std::string, int>::value_type("FULL", 2));

    return mapStringToInt.find(name)->second;
}
