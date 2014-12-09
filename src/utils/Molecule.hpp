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

#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include <iosfwd>
#include <vector>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Atom.hpp"
#include "Sphere.hpp"

enum rotorType {rtAsymmetric, rtSymmetric, rtSpherical, rtLinear, rtAtom};
const std::string rotorTypeList[] = {"Asymmetric", "Symmetric", "Spherical", "Linear", "Atom"};

/*! \file Molecule.hpp
 *  \class Molecule
 *  \brief Class representing a molecule or general aggregate of atoms.
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  This class is based on the similar class available in the Mints library
 *  of Psi4
 */

class Molecule {
private:
    /// The number of atoms in the molecule
    int nAtoms_;
    /// A vector of dimension (# atoms) containing the charges
    Eigen::VectorXd charges_;
    /// A vector of dimension (# atoms) containing the masses
    Eigen::VectorXd masses_;
    /// Molecular geometry, in cartesian coordinates. The dimensions are (# atoms * 3)
    /// Units are Bohr.
    Eigen::Matrix3Xd geometry_;
    /// A container for all the atoms composing the molecule
    std::vector<Atom> atoms_;
    /// A container for the spheres composing the molecule
    std::vector<Sphere> spheres_;
    /// The molecular rotor type
    rotorType rotor_;
public:
    Molecule() {}
    Molecule(int nat, const Eigen::VectorXd & chg, const Eigen::VectorXd & masses, 
             const Eigen::Matrix3Xd & geo, const std::vector<Atom> & at, const std::vector<Sphere> & sph); 
    Molecule(int nat, const std::vector<Sphere> & sph); 
    /// Copy constructor.
    Molecule(const Molecule &other);
    ~Molecule(){}
    
    int nAtoms() { return nAtoms_; }
    Eigen::VectorXd charges() { return charges_; }
    Eigen::VectorXd masses() { return masses_; }
    Eigen::Matrix3Xd geometry() { return geometry_; }
    std::vector<Atom> atoms() { return atoms_; }
    std::vector<Sphere> spheres() { return spheres_; }

    rotorType rotor();
    rotorType findRotorType();

    Eigen::Vector3d centerOfMass();
    Eigen::Matrix3d inertiaTensor();

    /*!
     * Given a vector, carries out translation of the molecule.
     * \param translationVector The translation vector.
     */
    void translate(const Eigen::Vector3d &translationVector);
    /// Performs translation to the Center of Mass Frame.
    void moveToCOM();

    /*!
     * Given a matrix, carries out rotation of the molecule.
     * \param rotationMatrix The matrix representing the rotation.
     */
    void rotate(const Eigen::Matrix3d &rotationMatrix);
    /// Performs rotation to the Principal Axes Frame.
    void moveToPAF();

    /// @{
    /// Operators
    /// Assignment operator.
    Molecule & operator=(const Molecule& other);
    friend std::ostream & operator<<(std::ostream &os, const Molecule &molecule);
    /// @}
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
};

#endif // MOLECULE_HPP
