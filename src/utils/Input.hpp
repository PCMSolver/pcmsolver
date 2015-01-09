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

#ifndef INPUT_HPP
#define INPUT_HPP

#include <vector>
#include <string>

#include "Config.hpp"

#include "Getkw.h"

#include "CavityData.hpp"
#include "GreenData.hpp"
#include "SolverData.hpp"
#include "InputManager.hpp"
#include "Molecule.hpp"
#include "Solvent.hpp"
#include "Sphere.hpp"

/*! \file Input.hpp
 *  \class Input
 *  \brief A wrapper class for the Getkw Library C++ bindings.
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  An Input object is to be used as the unique point of access to user-provided input:
 *   input ---> parsed input (Python script) ---> Input object (contains all the input data)
 *  Definition of input parameters is to be done in the Python script and in this class.
 *  They must be specified as private data members with public accessor methods (get-ters).
 *  Most of the data members are anyway accessed through the input wrapping struct-s
 *  In general, no mutator methods (set-ters) should be needed, exceptions to this rule
 *  should be carefully considered.
 */

class Input
{
public:
    /// Default constructor
    Input() {}
    /// Constructor from human-readable input file name
    Input(const std::string & filename);
    /// Constructor from host input structs
    Input(const cavityInput & cav, const solverInput & solv, const greenInput & green);
     
    /// Accessor methods
    
    /// Top-level section input
    std::string units() { return units_; }
    int CODATAyear() { return CODATAyear_; }
    /// @}

    /// Cavity section input
    std::string cavityType() { return type_; }
    bool scaling() { return scaling_; }
    std::string radiiSet() { return radiiSet_; }
    std::string mode() { return mode_; }
    std::vector<int> atoms() { return atoms_; }
    int atoms(int i) { return atoms_[i]; }
    std::vector<double> radii() { return radii_; }
    double radii(int i) { return radii_[i]; }
    std::vector<Sphere> spheres() { return spheres_; }
    Sphere spheres(int i) { return spheres_[i]; }
    Molecule molecule() { return molecule_; }
    /// This method sets the molecule and the list of spheres
    void molecule(const Molecule & m) { molecule_ = m; spheres_ = molecule_.spheres(); }
    /// @}
   
    /// Medium section input
    Solvent solvent() { return solvent_; }
    bool fromSolvent() { return hasSolvent_; }
    std::string solverType() { return solverType_; }
    int equationType() { return equationType_; }
    double correction() { return correction_; }
    bool hermitivitize() { return hermitivitize_; }
    /// @}

    /// Green's functin section input
    std::string greenInsideType() { return greenInsideType_; }
    std::string greenOutsideType() { return greenOutsideType_; }
    /// @}
    
    /// Keeps track of who did the parsing: the API or the host program
    std::string providedBy() { return providedBy_; }
    
    /// Get-ters for input wrapping structs
    cavityData cavityParams();
    greenData insideGreenParams();
    greenData outsideStaticGreenParams(); 
    greenData outsideDynamicGreenParams();
    /// @}

    /// Operators
    /// operator<<
    friend std::ostream & operator<<(std::ostream &os, const Input &input);
    /// @}
private:
    /*! Read Python-parsed input (API-side syntactic input parsing) into Input object
     */
    void reader(const char * pythonParsed);
    /*! Read host data structures (host-side syntactic input parsing) into Input object.
     *  It provides access to a **limited** number of options only, basically the ones
     *  that can be filled into the cavityInput, solverInput and greenInput data structures.
     *  Lengths and areas are **expected** to be in Angstrom/Angstrom^2 and will hence be converted
     *  to au/au^2.
     *  \note Specification of the solvent by name overrides any input given through the
     *  greenInput data structure!
     *  \warning The cavity can only be built in the "Implicit" mode, i.e. by grabbing the
     *  coordinates for the sphere centers from the host program.
     *  Atomic coordinates are **expected** to be in au!
     *  The "Atoms" and "Explicit" methods are only available using the explicit parsing
     *  by our Python script of a separate input file.
     */
    void reader(const cavityInput & cav, const solverInput & solv,
                const greenInput & green);
    /*! Perform semantic input parsing aka sanity check
     */
    void semanticCheck();

    /// Units of measure
    std::string units_;
    /// Year of the CODATA set to be used
    int CODATAyear_;
    /// The type of cavity
    std::string type_;
    /// Filename for the .npz cavity restart file
    std::string cavFilename_;
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
    /// Solvent probe radius
    double probeRadius_;
    /// Type of integrator for the diagonal of the boundary integral operators
    std::string integratorType_;
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
    /// Who performed the syntactic input parsing
    std::string providedBy_;
    /// Input wrapping struct for the cavity
    cavityData cavData_;
    /// Input wrapping struct for the Green's function inside
    greenData insideGreenData_;
    /// Input wrapping struct for the Green's function outside (static permittivity)
    greenData outsideStaticGreenData_;
    /// Input wrapping struct for the Green's function outside (dynamic permittivity)
    greenData outsideDynamicGreenData_;
};

/*!
 * A useful map to convert the Der string to an integer which will be passed to the Green's function CTOR.
 */
int derivativeTraits(const std::string & name);

/*!
 * A useful map to convert the EquationType string to an integer which will be passed to the Solver CTOR.
 */
int integralEquation(const std::string & name);

#endif // INPUT_HPP
