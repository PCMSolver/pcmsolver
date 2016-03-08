/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include <iosfwd>
#include <vector>
#include <string>

#include "Config.hpp"

#include <Eigen/Core>

class Element;

#include "Atom.hpp"
#include "Sphere.hpp"
#include "Symmetry.hpp"

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

class Molecule
{
private:
    /// The number of atoms in the molecule
    size_t nAtoms_;
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
    /// The molecular point group
    Symmetry pointGroup_;
public:
    /*! \brief Default constructor
     *  Initialize a dummy molecule, e.g. as placeholder, see Cavity.cpp loadCavity method
     */
    Molecule() { rotor_ = rtAsymmetric; pointGroup_ = buildGroup(0, 0, 0, 0); }
    /*! \brief Constructor from full molecular data
     *  \param[in] nat number of atoms
     *  \param[in] chg vector of atomic charges
     *  \param[in] masses vector of atomic masses
     *  \param[in] geo molecular geometry (format nat*3)
     *  \param[in] at vector of Atom objects
     *  \param[in] sph vector of Sphere objects
     *
     *  This initializes the molecule in C1 symmetry
     */
    Molecule(int nat, const Eigen::VectorXd & chg, const Eigen::VectorXd & masses,
             const Eigen::Matrix3Xd & geo, const std::vector<Atom> & at,
             const std::vector<Sphere> & sph);
    /*! \brief Constructor from full molecular data, plus number of generators and generators
     *  \param[in] nat number of atoms
     *  \param[in] chg vector of atomic charges
     *  \param[in] masses vector of atomic masses
     *  \param[in] geo molecular geometry (format nat*3)
     *  \param[in] at vector of Atom objects
     *  \param[in] sph vector of Sphere objects
     *  \param[in] nr_gen number of molecular point group generators
     *  \param[in] gen molecular point group generators
     *
     *  This initializes the molecule in the symmetry prescribed by nr_gen and gen.
     *  See documentation of the Symmetry object for the conventions.
     */
    Molecule(int nat, const Eigen::VectorXd & chg, const Eigen::VectorXd & masses,
             const Eigen::Matrix3Xd & geo, const std::vector<Atom> & at,
             const std::vector<Sphere> & sph,
             int nr_gen, int gen[3]);
    /*! \brief Constructor from full molecular data and point group
     *  \param[in] nat number of atoms
     *  \param[in] chg vector of atomic charges
     *  \param[in] masses vector of atomic masses
     *  \param[in] geo molecular geometry (format nat*3)
     *  \param[in] at vector of Atom objects
     *  \param[in] sph vector of Sphere objects
     *  \param[in] pg the molecular point group (a Symmetry object)
     *
     *  This initializes the molecule in the symmetry prescribed by pg.
     */
    Molecule(int nat, const Eigen::VectorXd & chg, const Eigen::VectorXd & masses,
             const Eigen::Matrix3Xd & geo, const std::vector<Atom> & at,
             const std::vector<Sphere> & sph,
             const Symmetry & pg);
    /*! \brief Constructor from list of spheres
     *  \param[in] sph  list of spheres
     *
     *  \warning This constructor is to be used **exclusively** when initializing the Molecule
     *  in EXPLICIT mode, i.e. when the user specifies explicitly spheres centers and radii.
     *
     *  Molecule is treated as an aggregate of spheres. We do not have information on the atomic
     *  species involved in the aggregate.
     *  Charges are set to 1.0; masses are set based on the radii; geometry is set from the list of spheres.
     *  All the atoms are dummy atoms. The point group is C1.
     */
    Molecule(const std::vector<Sphere> & sph);
    /// Copy constructor.
    Molecule(const Molecule &other);
    ~Molecule() {}

    size_t nAtoms() const { return nAtoms_; }
    Eigen::VectorXd charges() const { return charges_; }
    double charges(int i) const { return charges_(i); }
    Eigen::VectorXd masses() const { return masses_; }
    double masses(int i) const { return masses_(i); }
    Eigen::Matrix3Xd geometry() const { return geometry_; }
    double geometry(int i, int j) const { return geometry_(i, j); }
    std::vector<Atom> atoms() const { return atoms_; }
    Atom atoms(int i) const { return atoms_[i]; }
    std::vector<Sphere> spheres() const { return spheres_; }
    Sphere spheres(int i) const { return spheres_[i]; }

    rotorType rotor();
    rotorType findRotorType();

    Symmetry pointGroup() const { return pointGroup_; }
    void pointGroup(const Symmetry & pg) { pointGroup_ = pg; }

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

Eigen::VectorXd computeMEP(const Molecule & mol, const std::vector<Element> & el);
/*! \brief Compute MEP for a single point charge */
Eigen::VectorXd computeMEP(const std::vector<Element> & el, double charge = 1.0, const Eigen::Vector3d & origin = Eigen::Vector3d::Zero());

#endif // MOLECULE_HPP
