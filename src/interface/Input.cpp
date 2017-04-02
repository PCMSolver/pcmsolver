/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#include "Input.hpp"

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "Config.hpp"

#include <Eigen/Core>
#include "utils/getkw/Getkw.h"

#include <boost/algorithm/string.hpp>

#include "bi_operators/BIOperatorData.hpp"
#include "cavity/CavityData.hpp"
#include "green/GreenData.hpp"
#include "solver/SolverData.hpp"
#include "PCMInput.h"
#include "utils/Factory.hpp"
#include "utils/Solvent.hpp"
#include "utils/Sphere.hpp"

namespace pcm {

Input::Input(const std::string & filename) {
  reader(filename);
  semanticCheck();
}

Input::Input(const PCMInput & host_input) {
  reader(host_input);
  semanticCheck();
}

void Input::reader(const std::string & filename) {
  Getkw input_ = Getkw(filename, false, true);
  units_ = input_.getStr("UNITS");
  CODATAyear_ = input_.getInt("CODATA");
  initBohrToAngstrom(bohrToAngstrom, CODATAyear_);

  const Section & mol = input_.getSect("MOLECULE");
  MEPfromMolecule_ = true;
  if (mol.isDefined()) {
    geometry_ = mol.getDblVec("GEOMETRY");
    MEPfromMolecule_ = mol.getBool("MEP");
  }

  const Section & cavity = input_.getSect("CAVITY");

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
  radiiSet_ = boost::algorithm::to_upper_copy(cavity.getStr("RADIISET"));
  minimalRadius_ = cavity.getDbl("MINRADIUS");
  mode_ = boost::algorithm::to_upper_copy(cavity.getStr("MODE"));
  if (mode_ == "EXPLICIT") {
    std::vector<double> spheresInput = cavity.getDblVec("SPHERES");
    int j = 0;
    int nAtoms = int(spheresInput.size() / 4);
    for (int i = 0; i < nAtoms; ++i) {
      Eigen::Vector3d center;
      center = (Eigen::Vector3d() << spheresInput[j], spheresInput[j + 1], spheresInput[j + 2]).finished();
      Sphere sph(center, spheresInput[j + 3]);
      spheres_.push_back(sph);
      j += 4;
    }
    // Initialize molecule from spheres only when molecule section is absent
    if (!mol.isDefined())
      molecule_ = Molecule(spheres_);
  } else if (mode_ == "ATOMS") {
    atoms_ = cavity.getIntVec("ATOMS");
    radii_ = cavity.getDblVec("RADII");
  }

  // Get the contents of the Medium section
  const Section & medium = input_.getSect("MEDIUM");
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
    derivativeInsideType_ = detail::derivativeTraits(inside.getStr("DER"));
    epsilonInside_ = inside.getDbl("EPS");
    // Get the contents of the Green<outside> section...
    const Section & outside = medium.getSect("GREEN<OUTSIDE>");
    // ...and initialize the data members
    greenOutsideType_ = outside.getStr("TYPE");
    derivativeOutsideType_ = detail::derivativeTraits(outside.getStr("DER"));
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
    profileType_ = detail::profilePolicy(outside.getStr("PROFILE"));
    maxL_ = outside.getInt("MAXL");
  } else { // This part must be reviewed!! Some data members are not initialized...
    // Just initialize the solvent object in this class
    hasSolvent_ = true;
    if (solvents().find(name) == solvents().end()) {
      PCMSOLVER_ERROR("Solvent " + name + " NOT found!");
    } else {
      solvent_ = solvents()[name];
    }
    probeRadius_ = solvent_.probeRadius * angstromToBohr();
    // Specification of the solvent by name means isotropic PCM
    // We have to initialize the Green's functions data here, Solvent class
    // is an helper class and should not be used in the core classes.
    greenInsideType_ = "VACUUM";
    derivativeInsideType_ = detail::derivativeTraits("DERIVATIVE");
    epsilonInside_ = 1.0;
    greenOutsideType_ = "UNIFORMDIELECTRIC";
    derivativeOutsideType_ = detail::derivativeTraits("DERIVATIVE");
    epsilonStaticOutside_ = solvent_.epsStatic;
    epsilonDynamicOutside_ = solvent_.epsDynamic;
  }
  integratorType_ = medium.getStr("DIAGONALINTEGRATOR");
  integratorScaling_ = medium.getDbl("DIAGONALSCALING");

  solverType_ = medium.getStr("SOLVERTYPE");
  equationType_ = detail::integralEquation(medium.getStr("EQUATIONTYPE"));
  correction_ = medium.getDbl("CORRECTION");
  hermitivitize_ = medium.getBool("MATRIXSYMM");
  isDynamic_ = medium.getBool("NONEQUILIBRIUM");

  providedBy_ = std::string("API-side");
}

void Input::reader(const PCMInput & host_input) {
  CODATAyear_ = 2010;
  initBohrToAngstrom(bohrToAngstrom, CODATAyear_);

  type_ = detail::trim_and_upper(host_input.cavity_type);
  area_ = host_input.area * angstrom2ToBohr2();
  patchLevel_ = host_input.patch_level;
  coarsity_ = host_input.coarsity * angstromToBohr();
  minDistance_ = host_input.min_distance * angstromToBohr();
  derOrder_ = host_input.der_order;
  if (type_ == "RESTART") {
    cavFilename_ = detail::trim(host_input.restart_name); // No case conversion here!
  }

  scaling_ = host_input.scaling;
  radiiSet_ = detail::trim_and_upper(host_input.radii_set);
  if (radiiSet_ == "UFF") {
    radiiSetName_ = "UFF";
  } else if (radiiSet_ == "BONDI") {
    radiiSetName_ = "Bondi-Mantina";
  } else {
    radiiSetName_ = "Allinger's MM3";
  }
  minimalRadius_ = host_input.min_radius * angstromToBohr();
  mode_ = std::string("IMPLICIT");

  std::string name = detail::trim_and_upper(host_input.solvent);
  if (name.empty() || name == "EXPLICIT") {
    hasSolvent_ = false;
    // Get the probe radius
    probeRadius_ = host_input.probe_radius * angstromToBohr();
    // Get the contents of the Green<inside> section...
    // ...and initialize the data members
    greenInsideType_ = detail::trim_and_upper(host_input.inside_type);
    derivativeInsideType_ = detail::derivativeTraits("DERIVATIVE");
    epsilonInside_ = 1.0;
    // Get the contents of the Green<outside> section...
    // ...and initialize the data members
    greenOutsideType_ = detail::trim_and_upper(host_input.outside_type);
    derivativeOutsideType_ = detail::derivativeTraits("DERIVATIVE");
    epsilonStaticOutside_ = host_input.outside_epsilon;
    epsilonDynamicOutside_ = host_input.outside_epsilon;
    // Initialize interface parameters with bogus values
    epsilonStatic1_ = 0.0;
    epsilonDynamic1_ = 0.0;
    epsilonStatic2_ = 0.0;
    epsilonDynamic2_ = 0.0;
    center_ = 0.0;
    width_ = 0.0;
    origin_ = std::vector<double>(3, 0.0);
    profileType_ = 0;
    maxL_ = 0;
  } else { // This part must be reviewed!! Some data members are not initialized...
    // Just initialize the solvent object in this class
    hasSolvent_ = true;
    solvent_ = solvents()[name];
    probeRadius_ = solvent_.probeRadius * angstromToBohr();
    // Specification of the solvent by name means isotropic PCM
    // We have to initialize the Green's functions data here, Solvent class
    // is an helper class and should not be used in the core classes.
    greenInsideType_ = std::string("VACUUM");
    derivativeInsideType_ = detail::derivativeTraits("DERIVATIVE");
    epsilonInside_ = 1.0;
    greenOutsideType_ = std::string("UNIFORMDIELECTRIC");
    derivativeOutsideType_ = detail::derivativeTraits("DERIVATIVE");
    epsilonStaticOutside_ = solvent_.epsStatic;
    epsilonDynamicOutside_ = solvent_.epsDynamic;
  }

  integratorType_ = "COLLOCATION";
  integratorScaling_ = 1.07;

  solverType_ = detail::trim_and_upper(host_input.solver_type);
  std::string inteq = detail::trim_and_upper(host_input.equation_type);
  equationType_ = detail::integralEquation(inteq);
  correction_ = host_input.correction;
  hermitivitize_ = true;
  isDynamic_ = false;

  providedBy_ = std::string("host-side");
}

void Input::semanticCheck() {}

void Input::initMolecule() {
  // Gather information necessary to build molecule_
  // 1. number of atomic centers
  int nuclei = int(geometry_.size() / 4);
  // 2. position and charges of atomic centers
  Eigen::Matrix3Xd centers = Eigen::Matrix3Xd::Zero(3, nuclei);
  Eigen::VectorXd charges = Eigen::VectorXd::Zero(nuclei);
  int j = 0;
  for (int i = 0; i < nuclei; ++i) {
    centers.col(i) = (Eigen::Vector3d() << geometry_[j], geometry_[j + 1], geometry_[j + 2]).finished();
    charges(i) = geometry_[j + 3];
    j += 4;
  }
  // 3. list of atoms and list of spheres
  std::vector<Atom> radiiSet;
  std::vector<Atom> atoms;
  atoms.reserve(nuclei);
  // FIXME Code duplication in function initMolecule in interface/Meddle.cpp
  tie(radiiSetName_, radiiSet) = utils::bootstrapRadiiSet().create(radiiSet_);
  for (int i = 0; i < charges.size(); ++i) {
    int index = int(charges(i)) - 1;
    atoms.push_back(radiiSet[index]);
    if (scaling_)
      atoms[i].radiusScaling = 1.2;
  }
  // Based on the creation mode (Implicit or Atoms)
  // the spheres list might need postprocessing
  if (mode_ == "IMPLICIT" || mode_ == "ATOMS") {
    for (int i = 0; i < charges.size(); ++i) {
      // Convert to Bohr and multiply by scaling factor (alpha)
      double radius = atoms[i].radius * angstromToBohr() * atoms[i].radiusScaling;
      spheres_.push_back(Sphere(centers.col(i), radius));
    }
    if (mode_ == "ATOMS") {
      // Loop over the atomsInput array to get which atoms will have a user-given
      // radius
      for (size_t i = 0; i < atoms_.size(); ++i) {
        int index =
            atoms_[i] - 1; // -1 to go from human readable to machine readable
        // Put the new Sphere in place of the implicit-generated one
        spheres_[index] = Sphere(centers.col(index), radii_[i]);
      }
    }
  }

  // 4. masses
  Eigen::VectorXd masses = Eigen::VectorXd::Zero(nuclei);
  for (int i = 0; i < masses.size(); ++i) {
    masses(i) = atoms[i].mass;
  }
  // 5. molecular point group
  // FIXME currently hardcoded to C1

  // OK, now get molecule_
  molecule_ = Molecule(nuclei, charges, masses, centers, atoms, spheres_);
  // Check that all atoms have a radius attached
  std::vector<Atom>::const_iterator res =
      std::find_if(atoms.begin(), atoms.end(), invalid);
  if (res != atoms.end()) {
    std::cout << molecule_ << std::endl;
    PCMSOLVER_ERROR("Some atoms do not have a radius attached. Please specify a "
                    "radius for all atoms (see "
                    "http://pcmsolver.readthedocs.org/en/latest/users/input.html)!");
  }
}

CavityData Input::cavityParams() const {
  return CavityData(molecule_,
                    area_,
                    probeRadius_,
                    minDistance_,
                    derOrder_,
                    minimalRadius_,
                    patchLevel_,
                    coarsity_,
                    cavFilename_,
                    dyadicFilename_);
}

GreenData Input::insideGreenParams() const {
  return GreenData(derivativeInsideType_, profileType_, epsilonInside_);
}

GreenData Input::outsideStaticGreenParams() const {
  GreenData retval(derivativeOutsideType_, profileType_, epsilonStaticOutside_);
  if (not hasSolvent_) {
    retval.howProfile = profileType_;
    retval.epsilon1 = epsilonStatic1_;
    retval.epsilon2 = epsilonStatic2_;
    retval.center = center_;
    retval.width = width_;
    retval.origin = (Eigen::Vector3d() << origin_[0], origin_[1], origin_[2]).finished();
    retval.maxL = maxL_;
  }
  return retval;
}

GreenData Input::outsideDynamicGreenParams() const {
  GreenData retval(derivativeOutsideType_, profileType_, epsilonDynamicOutside_);
  if (not hasSolvent_) {
    retval.howProfile = profileType_;
    retval.epsilon1 = epsilonDynamic1_;
    retval.epsilon2 = epsilonDynamic2_;
    retval.center = center_;
    retval.width = width_;
    retval.origin = (Eigen::Vector3d() << origin_[0], origin_[1], origin_[2]).finished();
    retval.maxL = maxL_;
  }
  return retval;
}

SolverData Input::solverParams() const {
  return SolverData(correction_, equationType_, hermitivitize_);
}

BIOperatorData Input::integratorParams() const {
  return BIOperatorData(integratorScaling_);
}

namespace detail {
int derivativeTraits(const std::string & name) {
  static std::map<std::string, int> mapStringToInt;
  mapStringToInt.insert(std::map<std::string, int>::value_type("NUMERICAL", 0));
  mapStringToInt.insert(std::map<std::string, int>::value_type("DERIVATIVE", 1));
  mapStringToInt.insert(std::map<std::string, int>::value_type("GRADIENT", 2));
  mapStringToInt.insert(std::map<std::string, int>::value_type("HESSIAN", 3));

  return mapStringToInt.find(name)->second;
}

int profilePolicy(const std::string & name) {
  static std::map<std::string, int> mapStringToInt;
  mapStringToInt.insert(std::map<std::string, int>::value_type("TANH", 0));
  mapStringToInt.insert(std::map<std::string, int>::value_type("ERF", 1));

  return mapStringToInt.find(name)->second;
}

int integralEquation(const std::string & name) {
  static std::map<std::string, int> mapStringToInt;
  mapStringToInt.insert(std::map<std::string, int>::value_type("FIRSTKIND", 0));
  mapStringToInt.insert(std::map<std::string, int>::value_type("SECONDKIND", 1));
  mapStringToInt.insert(std::map<std::string, int>::value_type("FULL", 2));

  return mapStringToInt.find(name)->second;
}

std::string trim(const char * src) {
  std::string tmp(src);
  boost::algorithm::trim(tmp);
  return tmp;
}

std::string trim_and_upper(const char * src) {
  return boost::algorithm::to_upper_copy(trim(src));
}
} // namespace detail
} // namespace pcm
