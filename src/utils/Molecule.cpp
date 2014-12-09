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

#include "Molecule.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Atom.hpp"
#include "MathUtils.hpp"

Molecule::Molecule(int nat, const Eigen::VectorXd & chg, const Eigen::VectorXd & m, 
             const Eigen::MatrixX3d & geo, const std::vector<Atom> & at, const std::vector<Sphere> & sph)
	: nAtoms_(nat), charges_(chg), masses_(m), geometry_(geo), atoms_(at), spheres_(sph) 
{
    rotor_ = findRotorType();
}


Molecule::Molecule(const Molecule &other){
    *this = other;
}

Eigen::Vector3d Molecule::centerOfMass(){
    Eigen::Vector3d com;
    com << 0.0, 0.0, 0.0;
    for (int i = 0; i < nAtoms_; ++i){
        com += masses_(i) * atoms_[i].atomCoord();
    }
    com *= 1.0/masses_.sum();
    return com;
}

Eigen::Matrix3d Molecule::inertiaTensor(){
    Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();

    for (int i = 0; i < nAtoms_; ++i){
        // Diagonal
        inertia(0,0) += masses_(i) * (geometry_(i,1) * geometry_(i,1) + geometry_(i,2) * geometry_(i,2));
        inertia(1,1) += masses_(i) * (geometry_(i,0) * geometry_(i,0) + geometry_(i,2) * geometry_(i,2));
        inertia(2,2) += masses_(i) * (geometry_(i,0) * geometry_(i,0) + geometry_(i,1) * geometry_(i,1));

        // Off-diagonal
        inertia(0,1) -= masses_(i) * (geometry_(i,0) * geometry_(i,1));
        inertia(0,2) -= masses_(i) * (geometry_(i,0) * geometry_(i,2));
        inertia(1,2) -= masses_(i) * (geometry_(i,1) * geometry_(i,2));
    }
    // Now symmetrize
    hermitivitize(inertia);

    // Check elements for a numerical zero and make it a hard zero
    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            if (fabs(inertia(i,j)) < 1.0e-14) {
                inertia(i,j) = 0.0;
            }
        }
    }

    return inertia;
}

rotorType Molecule::findRotorType(){
    rotorType type;
    if (nAtoms_ == 1) {
        type = rtAtom;
    } else {
        // Get inertia tensor
        Eigen::Matrix3d inertia = inertiaTensor();
        // Diagonalize inertia tensor V^t * I * V
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(inertia);
        if (eigenSolver.info() != Eigen::Success) abort();
        // Determine the degeneracy of the eigenvalues.
        int deg = 0;
        double tmp, abs, rel;
        for (int i = 0; i < 2; ++i){
            for (int j = i + 1; j < 3 && deg < 2; ++j){ // Check i and j != i
                abs = fabs(eigenSolver.eigenvalues()[i] - eigenSolver.eigenvalues()[j]);
                tmp = eigenSolver.eigenvalues()[j]; // Because the eigenvalues are already in ascending order.
                if (abs > 1.0e-14) {
                    rel = abs/tmp;
                } else {
                    rel = 0.0;
                }
                if (rel < 1.0e-8) {
                    ++deg;
                }
            }
        }
        // Get the rotor type based on the degeneracy.
        if (eigenSolver.eigenvalues()[0] == 0.0) {
            type = rtLinear;
        } else if (deg == 2) {
            type = rtSpherical;
        } else if (deg == 1) { // We do not distinguish between prolate and oblate.
            type = rtSymmetric;
        } else {
            type = rtAsymmetric;
        }
    }

    return type;
}

void Molecule::translate(const Eigen::Vector3d &translationVector){
    // Translate the geometry_ matrix and update the geometric data in atoms_.
    for (int i = 0; i < nAtoms_; ++i){
        geometry_.row(i) -= translationVector;
	Eigen::Vector3d tmp = geometry_.row(i).transpose();
        atoms_[i].atomCoord(tmp);
    }
}

void Molecule::moveToCOM(){
    Eigen::Vector3d com = centerOfMass();
    this->translate(com);
}

void Molecule::rotate(const Eigen::Matrix3d &rotationMatrix){
    // Rotate the geometry_ matrix and update the geometric data in atoms_.
    geometry_ *= rotationMatrix; // The power of Eigen: geometry_ = geometry_ * rotationMatrix;
    for (int i = 0; i < nAtoms_; ++i){
	Eigen::Vector3d tmp = geometry_.row(i).transpose();
        atoms_[i].atomCoord(tmp);
    }
}

void Molecule::moveToPAF(){
    Eigen::Matrix3d inertia = inertiaTensor();
    // Diagonalize inertia tensor V^t * I * V
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(inertia);
    if (eigenSolver.info() != Eigen::Success) abort();
    // Rotate to Principal Axes Frame
    this->rotate(eigenSolver.eigenvectors());
    std::cout << eigenSolver.eigenvalues() << std::endl;
}

Molecule& Molecule::operator=(const Molecule& other){
    // Self assignment is bad
    if (this == &other)
        return *this;

    nAtoms_ = other.nAtoms_;
    charges_ = other.charges_;
    masses_ = other.masses_;
    geometry_ = other.geometry_;
    atoms_ = other.atoms_;
    spheres_ = other.spheres_;
    rotor_ = other.rotor_;

    return *this;
}

std::ostream & operator<<(std::ostream &os, const Molecule &m){
    // Declare formatting of Eigen output.
    std::string sep = "                  ";
    Eigen::IOFormat CleanFmt(Eigen::FullPrecision, Eigen::DontAlignCols, sep, "\n", "", "");

    os << "Rotor type: " << rotorTypeList[m.rotor_] << std::endl;
    if (m.nAtoms_ != 0) {
        os << "       Center              X                  Y                   Z       " << std::endl;
        os << "    ------------   -----------------  -----------------  -----------------" << std::endl;
        for (int i = 0; i < m.nAtoms_; ++i){
            os << std::setw(10) << m.atoms_[i].atomSymbol() << std::setw(15) <<m.geometry_.row(i).format(CleanFmt) << std::endl;
        }
    } else {
        os << "  No atoms in this molecule!" << std::endl;
    }

    return os;
}
