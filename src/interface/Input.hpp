/*
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

#pragma once

#include <string>
#include <vector>

#include "Config.hpp"

#include "utils/getkw/Getkw.h"

/*! \file Input.hpp */

/*! \struct PCMInput
 */

struct PCMInput;

namespace pcm {
struct BIOperatorData;
struct CavityData;
struct GreenData;
struct SolverData;
} // namespace pcm

#include "utils/ChargeDistribution.hpp"
#include "utils/Molecule.hpp"
#include "utils/Solvent.hpp"
#include "utils/Sphere.hpp"

namespace pcm {
using utils::ChargeDistribution;
using utils::Solvent;
using utils::Sphere;

/*! \class Input
 *  \brief A wrapper class for the Getkw Library C++ bindings.
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  An Input object is to be used as the unique point of access to user-provided
 *input:
 *   input ---> parsed input (Python script) ---> Input object (contains all the
 *input data)
 *  Definition of input parameters is to be done in the Python script and in this
 *class.
 *  They must be specified as private data members with public accessor methods
 *(get-ters).
 *  Most of the data members are anyway accessed through the input wrapping struct-s
 *  In general, no mutator methods (set-ters) should be needed, exceptions to this
 *rule
 *  should be carefully considered.
 */
class Input {
public:
  /// Default constructor
  Input() {}
  /// Constructor from human-readable input file name
  Input(const std::string & filename);
  /// Constructor from host input structs
  Input(const PCMInput & host_input);

  /// Accessor methods

  /// Top-level section input
  std::string units() const { return units_; }
  int CODATAyear() const { return CODATAyear_; }
  /// @}

  /// Cavity section input
  std::string cavityType() const { return type_; }
  bool scaling() const { return scaling_; }
  std::string radiiSet() const { return radiiSet_; }
  std::string radiiSetName() const { return radiiSetName_; }
  std::string mode() const { return mode_; }
  std::vector<int> atoms() const { return atoms_; }
  int atoms(size_t i) const { return atoms_[i]; }
  std::vector<double> radii() const { return radii_; }
  double radii(size_t i) const { return radii_[i]; }
  std::vector<Sphere> spheres() const { return spheres_; }
  Sphere spheres(int i) const { return spheres_[i]; }
  Molecule molecule() const { return molecule_; }
  /// This method sets the molecule and the list of spheres
  void molecule(const Molecule & m) {
    molecule_ = m;
    spheres_ = molecule_.spheres();
  }
  void initMolecule();
  /// @}

  /// Medium section input
  Solvent solvent() const { return solvent_; }
  bool fromSolvent() const { return hasSolvent_; }
  std::string solverType() const { return solverType_; }
  int equationType() const { return equationType_; }
  double correction() const { return correction_; }
  bool hermitivitize() const { return hermitivitize_; }
  bool isDynamic() const { return isDynamic_; }
  std::string integratorType() const { return integratorType_; }
  double integratorScaling() const { return integratorScaling_; }
  /// @}

  /// Green's function section input
  std::string greenInsideType() const { return greenInsideType_; }
  std::string greenOutsideType() const { return greenOutsideType_; }
  /// @}

  /// Keeps track of who did the parsing: the API or the host program
  std::string providedBy() const { return providedBy_; }

  /// Get-ters for input wrapping structs
  CavityData cavityParams() const;
  GreenData insideGreenParams() const;
  GreenData outsideStaticGreenParams() const;
  GreenData outsideDynamicGreenParams() const;
  SolverData solverParams() const;
  BIOperatorData integratorParams() const;
  /// @}

  ChargeDistribution multipoles() const { return multipoles_; }
  bool MEPfromMolecule() { return MEPfromMolecule_; }

  /// Operators
  /// operator<<
  friend std::ostream & operator<<(std::ostream & os, const Input & input);
  /// @}
private:
  void reader(const std::string & filename);
  /*! Read host data structures (host-side syntactic input parsing) into Input
   * object.
   *  It provides access to a **limited** number of options only, basically the ones
   *  that can be filled into the cavityInput, solverInput and greenInput data
   * structures.
   *  Lengths and areas are **expected** to be in Angstrom/Angstrom^2 and will hence
   * be converted
   *  to au/au^2.
   *  \note Specification of the solvent by name overrides any input given through
   * the
   *  greenInput data structure!
   *  \warning The cavity can only be built in the "Implicit" mode, i.e. by grabbing
   * the
   *  coordinates for the sphere centers from the host program.
   *  Atomic coordinates are **expected** to be in au!
   *  The "Atoms" and "Explicit" methods are only available using the explicit
   * parsing
   *  by our Python script of a separate input file.
   */
  void reader(const PCMInput & host_input);
  /*! Perform semantic input parsing aka sanity check */
  void semanticCheck() attribute(const);

  /// Units of measure
  std::string units_;
  /// Year of the CODATA set to be used
  int CODATAyear_;
  /// The type of cavity
  std::string type_;
  /// Filename for the .npz cavity restart file
  std::string cavFilename_;
  /// Filename for the wavelet cavity dyadic file
  std::string dyadicFilename_;
  /// Wavelet cavity patch level
  int patchLevel_;
  /// Wavelet cavity coarsity
  double coarsity_;
  /// GePol and TsLess cavities average element area
  double area_;
  /// TsLess cavity minimal distance between sampling points
  double minDistance_;
  /// TsLess cavity maximum derivative order of switch function
  int derOrder_;
  /// Whether the radii should be scaled by 1.2
  bool scaling_;
  /// The set of radii to be used
  std::string radiiSet_;
  /// Collects info on atomic radii set
  std::string radiiSetName_;
  /// Minimal radius of an added sphere
  double minimalRadius_;
  /// How the API should get the coordinates of the sphere centers
  std::string mode_;
  /// List of selected atoms with custom radii
  std::vector<int> atoms_;
  /// List of radii attached to the selected atoms
  std::vector<double> radii_;
  /// List of spheres for fully custom cavity generation
  std::vector<Sphere> spheres_;
  /// Molecule or atomic aggregate
  Molecule molecule_;
  /// The solvent for a vacuum/uniform dielectric run
  Solvent solvent_;
  /// Whether the medium was initialized from a solvent object
  bool hasSolvent_;
  /// The solver type
  std::string solverType_;
  /// The integral equation type (wavelet solvers)
  int equationType_;
  /// Correction factor (C-PCM)
  double correction_;
  /// Whether the PCM matrix should be hermitivitized (collocation solvers)
  bool hermitivitize_;
  /// Whether the dynamic PCM matrix should be used
  bool isDynamic_;
  /// Solvent probe radius
  double probeRadius_;
  /// Type of integrator for the diagonal of the boundary integral operators
  std::string integratorType_;
  /// Scaling factor for the diagonal of the approximate collocation boundary
  /// integral operators
  double integratorScaling_;
  /// The Green's function type inside the cavity
  std::string greenInsideType_;
  /// The Green's function type outside the cavity
  std::string greenOutsideType_;
  /// How to calculate Green's function derivatives inside the cavity
  int derivativeInsideType_;
  /// How to calculate Green's function derivatives outside the cavity
  int derivativeOutsideType_;
  /// Permittivity inside the cavity
  double epsilonInside_;
  /// Static permittivity outside the cavity
  double epsilonStaticOutside_;
  /// Dynamic permittivity outside the cavity
  double epsilonDynamicOutside_;
  /// Real part of the metal NP permittivity
  double epsilonReal_;
  /// Imaginary part of the metal NP permittivity
  double epsilonImaginary_;
  /// Center of the spherical metal NP
  std::vector<double> spherePosition_;
  /// Radius of the spherical metal NP
  double sphereRadius_;
  /// Diffuse interface: static permittivity inside the interface
  double epsilonStatic1_;
  /// Diffuse interface: dynamic permittivity inside the interface
  double epsilonDynamic1_;
  /// Diffuse interface: static permittivity outside the interface
  double epsilonStatic2_;
  /// Diffuse interface: dynamic permittivity outside the interface
  double epsilonDynamic2_;
  /// Center of the diffuse interface
  double center_;
  /// Width of the diffuse interface
  double width_;
  /// Profile chosen for the diffuse interface
  int profileType_;
  /// Maximum angular momentum
  int maxL_;
  /// Center of the dielectric sphere
  std::vector<double> origin_;
  /// Molecular geometry
  std::vector<double> geometry_;
  /// Whether to calculate the MEP from the molecular geometry
  bool MEPfromMolecule_;
  /// Classical charge distribution of point multipoles
  ChargeDistribution multipoles_;
  /// Who performed the syntactic input parsing
  std::string providedBy_;
};

namespace detail {
/*! A useful map to convert the Der string to an integer which will be passed to the
 * Green's function CTOR. */
int derivativeTraits(const std::string & name);

/*! A useful map to convert from a string specifying the dielectric profile to an
 * integer
 *  which will be passed to the Green's function CTOR.
 */
int profilePolicy(const std::string & name);

/*! A useful map to convert the EquationType string to an integer which will be
 * passed to the Solver CTOR. */
int integralEquation(const std::string & name);

std::string trim(const char * src);

std::string trim_and_upper(const char * src);
} // namespace detail
} // namespace pcm
