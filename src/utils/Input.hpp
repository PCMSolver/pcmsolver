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

#include "Solvent.hpp"
#include "Sphere.hpp"

/*! \file Input.hpp
 *  \class Input
 *  \brief A wrapper class for the Getkw Library C++ bindings.
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  Implemented as a Singleton: only one instance of it is needed and therefore admitted.
 *  Constructors (default, non-default and copy), destructor and assignment operator are
 *  made private. We use the lazy initialization idiom to initialize the unique Input
 *  object \cite Alexandrescu2001
 *  An Input object is to be used as the unique point of access to user-provided input:
 *   input ---> parsed input (Python script) ---> Input object (contains all the input data)
 *  Definition of input parameters is to be done in the Python script and in this class.
 *  At runtime a unique Input object in created, this object contains all input parameters
 *  processed as to be directly usable in the module.
 *  They must be specified as private data members with public accessor methods (get-ters).
 *  In general, no mutator methods (set-ters) should be needed, exceptions to this rule
 *  should be carefully considered.
 */

class Input
{
private:
    Input();
    Input(const Input &other);
    Input& operator=(const Input &other);
    ~Input() {}
    int CODATAyear_;
    std::string cavFilename;
    std::string type;
    int patchLevel;
    double coarsity;
    double area;
    double minDistance;
    int derOrder;
    bool scaling;
    std::string radiiSet;
    double minimalRadius;
    std::string mode;
    std::vector<int> atoms;
    std::vector<double> radii;
    std::vector<Sphere> spheres;
    Solvent solvent;
    bool hasSolvent;
    std::string solverType;
    int equationType;
    double correction;
    bool hermitivitize_;
    double probeRadius;
    std::string greenInsideType;
    std::string greenOutsideType;
    int derivativeInsideType;
    int derivativeOutsideType;
    double epsilonInside;
    double epsilonOutside;
    double epsilonReal;
    double epsilonImaginary;
    std::vector<double> spherePosition;
    double sphereRadius;
public:
    static Input& TheInput() {
        static Input obj;
        return obj;
    }
    // Accessor methods
    int CODATAyear() { return CODATAyear_; } 
    // Cavity section input
    std::string getCavityFilename() { return cavFilename; }
    std::string getCavityType() { return type; }
    int getPatchLevel() { return patchLevel; }
    double getCoarsity() { return coarsity; }
    double getArea() { return area; }
    double getMinDistance() { return minDistance; }
    int getDerOrder() { return derOrder; }
    bool getScaling() { return scaling; }
    std::string getRadiiSet() { return radiiSet; }
    double getMinimalRadius() { return minimalRadius; }
    std::string getMode() { return mode; }
    std::vector<int> getAtoms() { return atoms; }
    std::vector<double> getRadii() { return radii; }
    std::vector<Sphere> getSpheres() { return spheres; }
    void setSpheres(const std::vector<Sphere> & _spheres) { spheres = _spheres; }
    // Medium section input
    Solvent getSolvent() { return solvent; }
    bool fromSolvent() { return hasSolvent; }
    std::string getSolverType() { return solverType; }
    int getEquationType() { return equationType; }
    double getCorrection() { return correction; }
    bool hermitivitize() { return hermitivitize_; }
    double getProbeRadius() { return probeRadius; }
    std::string getGreenInsideType() { return greenInsideType; }
    std::string getGreenOutsideType() { return greenOutsideType; }
    int getDerivativeInsideType() { return derivativeInsideType; }
    int getDerivativeOutsideType() { return derivativeOutsideType; }
    double getEpsilonInside() { return epsilonInside; }
    double getEpsilonOutside() { return epsilonOutside; }
    double getEpsilonReal() { return epsilonReal; }
    double getEpsilonImaginary() { return epsilonImaginary; }
    std::vector<double> getSpherePosition() { return spherePosition; }
    double getSphereRadius() { return sphereRadius; }
    /// Operators
    /// operator<<
    friend std::ostream & operator<<(std::ostream &os, const Input &input);
    /// @}
};

#endif // INPUT_HPP
