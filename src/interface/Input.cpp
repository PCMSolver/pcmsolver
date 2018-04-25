/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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

#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Config.hpp"

#include "utils/getkw/Getkw.h"
#include <Eigen/Core>

#include "PCMInput.h"
#include "bi_operators/BIOperatorData.hpp"
#include "cavity/CavityData.hpp"
#include "green/GreenData.hpp"
#include "solver/SolverData.hpp"
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

  cavityType_ = cavity.getStr("TYPE");
  area_ = cavity.getDbl("AREA");
  if (cavityType_ == "RESTART") {
    cavFilename_ = cavity.getStr("NPZFILE");
  }

  scaling_ = cavity.getBool("SCALING");
  radiiSet_ = detail::uppercase(cavity.getStr("RADIISET"));
  minimalRadius_ = cavity.getDbl("MINRADIUS");
  mode_ = detail::uppercase(cavity.getStr("MODE"));
  if (mode_ == "EXPLICIT") {
    std::vector<double> spheresInput = cavity.getDblVec("SPHERES");
    int j = 0;
    int nAtoms = int(spheresInput.size() / 4);
    for (int i = 0; i < nAtoms; ++i) {
      Eigen::Vector3d center;
      center = (Eigen::Vector3d() << spheresInput[j],
                spheresInput[j + 1],
                spheresInput[j + 2])
                   .finished();
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
    greenInsideType_ = inside.getStr("TYPE") + "_" + inside.getStr("DER");
    epsilonInside_ = inside.getDbl("EPS");
    // Get the contents of the Green<outside> section...
    const Section & outside = medium.getSect("GREEN<OUTSIDE>");
    // ...and initialize the data members
    greenOutsideType_ = outside.getStr("TYPE") + "_" + outside.getStr("DER");
    epsilonStaticOutside_ = outside.getDbl("EPS");
    epsilonDynamicOutside_ = outside.getDbl("EPSDYN");
    epsilonStatic1_ = outside.getDbl("EPS1");
    epsilonDynamic1_ = outside.getDbl("EPSDYN1");
    epsilonStatic2_ = outside.getDbl("EPS2");
    epsilonDynamic2_ = outside.getDbl("EPSDYN2");
    center_ = outside.getDbl("CENTER");
    width_ = outside.getDbl("WIDTH");
    origin_ = outside.getDblVec("INTERFACEORIGIN");
    if (outside.getStr("TYPE") == "SPHERICALDIFFUSE") {
      greenOutsideType_ += "_" + outside.getStr("PROFILE");
    }
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
    greenInsideType_ = "VACUUM_DERIVATIVE";
    epsilonInside_ = 1.0;
    greenOutsideType_ = "UNIFORMDIELECTRIC_DERIVATIVE";
    epsilonStaticOutside_ = solvent_.epsStatic;
    epsilonDynamicOutside_ = solvent_.epsDynamic;
  }
  integratorType_ = medium.getStr("DIAGONALINTEGRATOR");
  integratorScaling_ = medium.getDbl("DIAGONALSCALING");

  solverType_ = medium.getStr("SOLVERTYPE");
  correction_ = medium.getDbl("CORRECTION");
  hermitivitize_ = medium.getBool("MATRIXSYMM");
  isDynamic_ = medium.getBool("NONEQUILIBRIUM");

  const Section & chgdist = input_.getSect("CHARGEDISTRIBUTION");
  if (chgdist.isDefined()) {
    // Set monopoles
    if (chgdist.getKey<std::vector<double> >("MONOPOLES").isDefined()) {
      std::vector<double> mono = chgdist.getDblVec("MONOPOLES");
      int j = 0;
      int n = int(mono.size() / 4);
      multipoles_.monopoles = Eigen::VectorXd::Zero(n);
      multipoles_.monopolesSites = Eigen::Matrix3Xd::Zero(3, n);
      for (int i = 0; i < n; ++i) {
        multipoles_.monopolesSites.col(i) =
            (Eigen::Vector3d() << mono[j], mono[j + 1], mono[j + 2]).finished();
        multipoles_.monopoles(i) = mono[j + 3];
        j += 4;
      }
    }
    // Set dipoles
    if (chgdist.getKey<std::vector<double> >("DIPOLES").isDefined()) {
      std::vector<double> dipo = chgdist.getDblVec("DIPOLES");
      int j = 0;
      int n = int(dipo.size() / 6);
      multipoles_.dipoles = Eigen::Matrix3Xd::Zero(3, n);
      multipoles_.dipolesSites = Eigen::Matrix3Xd::Zero(3, n);
      for (int i = 0; i < n; ++i) {
        multipoles_.dipolesSites.col(i) =
            (Eigen::Vector3d() << dipo[j], dipo[j + 1], dipo[j + 2]).finished();
        multipoles_.dipoles.col(i) =
            (Eigen::Vector3d() << dipo[j + 3], dipo[j + 4], dipo[j + 5]).finished();
        j += 6;
      }
    }
  }

  providedBy_ = std::string("API-side");
}

void Input::reader(const PCMInput & host_input) {
  CODATAyear_ = 2010;
  initBohrToAngstrom(bohrToAngstrom, CODATAyear_);

  cavityType_ = detail::trim_and_upper(host_input.cavity_type);
  area_ = host_input.area * angstrom2ToBohr2();
  if (cavityType_ == "RESTART") {
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
    greenInsideType_ =
        detail::trim_and_upper(host_input.inside_type) + "_DERIVATIVE";
    epsilonInside_ = 1.0;
    // Get the contents of the Green<outside> section...
    // ...and initialize the data members
    greenOutsideType_ =
        detail::trim_and_upper(host_input.outside_type) + "_DERIVATIVE";
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
    maxL_ = 0;
  } else { // This part must be reviewed!! Some data members are not initialized...
    // Just initialize the solvent object in this class
    hasSolvent_ = true;
    solvent_ = solvents()[name];
    probeRadius_ = solvent_.probeRadius * angstromToBohr();
    // Specification of the solvent by name means isotropic PCM
    // We have to initialize the Green's functions data here, Solvent class
    // is an helper class and should not be used in the core classes.
    greenInsideType_ = std::string("VACUUM_DERIVATIVE");
    epsilonInside_ = 1.0;
    greenOutsideType_ = std::string("UNIFORMDIELECTRIC_DERIVATIVE");
    epsilonStaticOutside_ = solvent_.epsStatic;
    epsilonDynamicOutside_ = solvent_.epsDynamic;
  }

  integratorType_ = "COLLOCATION";
  integratorScaling_ = 1.07;

  solverType_ = detail::trim_and_upper(host_input.solver_type);
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
    centers.col(i) =
        (Eigen::Vector3d() << geometry_[j], geometry_[j + 1], geometry_[j + 2])
            .finished();
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
  return CavityData(
      cavityType_, molecule_, area_, probeRadius_, minimalRadius_, cavFilename_);
}

GreenData Input::insideGreenParams() const {
  return GreenData(greenInsideType_, epsilonInside_);
}

GreenData Input::outsideStaticGreenParams() const {
  GreenData retval(greenOutsideType_, epsilonStaticOutside_);
  if (not hasSolvent_) {
    retval.epsilon1 = epsilonStatic1_;
    retval.epsilon2 = epsilonStatic2_;
    retval.center = center_;
    retval.width = width_;
    retval.origin =
        (Eigen::Vector3d() << origin_[0], origin_[1], origin_[2]).finished();
    retval.maxL = maxL_;
  }
  return retval;
}

GreenData Input::outsideDynamicGreenParams() const {
  GreenData retval(greenOutsideType_, epsilonDynamicOutside_);
  if (not hasSolvent_) {
    retval.epsilon1 = epsilonDynamic1_;
    retval.epsilon2 = epsilonDynamic2_;
    retval.center = center_;
    retval.width = width_;
    retval.origin =
        (Eigen::Vector3d() << origin_[0], origin_[1], origin_[2]).finished();
    retval.maxL = maxL_;
  }
  return retval;
}

SolverData Input::solverParams() const {
  return SolverData(solverType_, correction_, hermitivitize_);
}

BIOperatorData Input::integratorParams() const {
  return BIOperatorData(integratorType_, integratorScaling_);
}

namespace detail {
#ifndef HAS_CXX11
std::string left_trim(std::string s) {
  s.erase(s.begin(),
          std::find_if(
              s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}

std::string right_trim(std::string s) {
  s.erase(std::find_if(
              s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace)))
              .base(),
          s.end());
  return s;
}
#else  /* HAS_CXX11 */
std::string left_trim(std::string s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
            return !std::isspace(ch);
          }));
  return s;
}

std::string right_trim(std::string s) {
  s.erase(
      std::find_if(s.rbegin(), s.rend(), [](int ch) { return !std::isspace(ch); })
          .base(),
      s.end());
  return s;
}
#endif /* HAS_CXX11 */

std::string left_trim(const char * src) {
  std::string tmp(src);
  return left_trim(tmp);
}

std::string right_trim(const char * src) {
  std::string tmp(src);
  return right_trim(tmp);
}

std::string trim(std::string s) {
  left_trim(s);
  right_trim(s);
  return s;
}

std::string trim(const char * src) {
  std::string tmp(src);
  return trim(tmp);
}

std::string uppercase(std::string s) {
  std::transform(
      s.begin(), s.end(), s.begin(), std::ptr_fun<int, int>(std::toupper));
  return s;
}

std::string uppercase(const char * src) {
  std::string tmp(src);
  return uppercase(tmp);
}

std::string trim_and_upper(const char * src) {
  std::string tmp(src);
  return uppercase(trim(tmp));
}
} // namespace detail
} // namespace pcm
