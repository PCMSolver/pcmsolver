/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
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

#include "Molecule.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "Atom.hpp"
#include "MathUtils.hpp"
#include "Symmetry.hpp"
#include "cavity/Element.hpp"

namespace pcm {
using cavity::Element;
using utils::hermitivitize;

Molecule::Molecule(int nat,
                   const Eigen::VectorXd & chg,
                   const Eigen::VectorXd & m,
                   const Eigen::Matrix3Xd & geo,
                   const std::vector<Atom> & at,
                   const std::vector<Sphere> & sph)
    : nAtoms_(nat),
      charges_(chg),
      masses_(m),
      geometry_(geo),
      atoms_(at),
      spheres_(sph) {
  rotor_ = findRotorType();
  pointGroup_ = buildGroup(0, 0, 0, 0);
}

Molecule::Molecule(int nat,
                   const Eigen::VectorXd & chg,
                   const Eigen::VectorXd & m,
                   const Eigen::Matrix3Xd & geo,
                   const std::vector<Atom> & at,
                   const std::vector<Sphere> & sph,
                   int nr_gen,
                   int gen[3])
    : nAtoms_(nat),
      charges_(chg),
      masses_(m),
      geometry_(geo),
      atoms_(at),
      spheres_(sph) {
  rotor_ = findRotorType();
  pointGroup_ = buildGroup(nr_gen, gen[0], gen[1], gen[2]);
}

Molecule::Molecule(int nat,
                   const Eigen::VectorXd & chg,
                   const Eigen::VectorXd & m,
                   const Eigen::Matrix3Xd & geo,
                   const std::vector<Atom> & at,
                   const std::vector<Sphere> & sph,
                   const Symmetry & pg)
    : nAtoms_(nat),
      charges_(chg),
      masses_(m),
      geometry_(geo),
      atoms_(at),
      spheres_(sph),
      pointGroup_(pg) {
  rotor_ = findRotorType();
}

Molecule::Molecule(const std::vector<Sphere> & sph)
    : nAtoms_(sph.size()), spheres_(sph) {
  charges_ = Eigen::VectorXd::Ones(nAtoms_);
  masses_.resize(nAtoms_);
  geometry_.resize(Eigen::NoChange, nAtoms_);
  for (size_t i = 0; i < nAtoms_; ++i) {
    masses_(i) = spheres_[i].radius;
    geometry_.col(i) = spheres_[i].center;
    double charge = charges_(i);
    double mass = masses_(i);
    atoms_.push_back(Atom("Dummy", "Du", charge, mass, mass, geometry_.col(i)));
  }
  rotor_ = findRotorType();
  pointGroup_ = buildGroup(0, 0, 0, 0);
}

Molecule::Molecule(const Molecule & other) { *this = other; }

Eigen::Vector3d Molecule::centerOfMass() {
  Eigen::Vector3d com;
  com << 0.0, 0.0, 0.0;
  for (size_t i = 0; i < nAtoms_; ++i) {
    com += masses_(i) * atoms_[i].position;
  }
  com *= 1.0 / masses_.sum();
  return com;
}

Eigen::Matrix3d Molecule::inertiaTensor() {
  Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();

  for (size_t i = 0; i < nAtoms_; ++i) {
    // Diagonal
    inertia(0, 0) += masses_(i) * (geometry_(1, i) * geometry_(1, i) +
                                   geometry_(2, i) * geometry_(2, i));
    inertia(1, 1) += masses_(i) * (geometry_(0, i) * geometry_(0, i) +
                                   geometry_(2, i) * geometry_(2, i));
    inertia(2, 2) += masses_(i) * (geometry_(0, i) * geometry_(0, i) +
                                   geometry_(1, i) * geometry_(1, i));

    // Off-diagonal
    inertia(0, 1) -= masses_(i) * (geometry_(0, i) * geometry_(1, i));
    inertia(0, 2) -= masses_(i) * (geometry_(0, i) * geometry_(2, i));
    inertia(1, 2) -= masses_(i) * (geometry_(1, i) * geometry_(2, i));
  }
  // Now symmetrize
  hermitivitize(inertia);

  // Check elements for a numerical zero and make it a hard zero
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (fabs(inertia(i, j)) < 1.0e-14) {
        inertia(i, j) = 0.0;
      }
    }
  }

  return inertia;
}

rotorType Molecule::findRotorType() {
  rotorType type;
  if (nAtoms_ == 1) {
    type = rtAtom;
  } else {
    // Get inertia tensor
    Eigen::Matrix3d inertia = inertiaTensor();
    // Diagonalize inertia tensor V^t * I * V
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(inertia);
    if (eigenSolver.info() != Eigen::Success)
      abort();
    // Determine the degeneracy of the eigenvalues.
    int deg = 0;
    double tmp, abs, rel;
    for (int i = 0; i < 2; ++i) {
      for (int j = i + 1; j < 3 && deg < 2; ++j) { // Check i and j != i
        abs = fabs(eigenSolver.eigenvalues()[i] - eigenSolver.eigenvalues()[j]);
        tmp = eigenSolver.eigenvalues()[j]; // Because the eigenvalues are already in
                                            // ascending order.
        if (abs > 1.0e-14) {
          rel = abs / tmp;
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

void Molecule::translate(const Eigen::Vector3d & translationVector) {
  // Translate the geometry_ matrix and update the geometric data in atoms_.
  for (size_t i = 0; i < nAtoms_; ++i) {
    geometry_.col(i) -= translationVector;
    Eigen::Vector3d tmp = geometry_.col(i);
    atoms_[i].position = tmp;
  }
}

void Molecule::moveToCOM() {
  Eigen::Vector3d com = centerOfMass();
  this->translate(com);
}

void Molecule::rotate(const Eigen::Matrix3d & rotationMatrix) {
  // Rotate the geometry_ matrix and update the geometric data in atoms_.
  geometry_ *=
      rotationMatrix; // The power of Eigen: geometry_ = geometry_ * rotationMatrix;
  for (size_t i = 0; i < nAtoms_; ++i) {
    Eigen::Vector3d tmp = geometry_.col(i);
    atoms_[i].position = tmp;
  }
}

void Molecule::moveToPAF() {
  Eigen::Matrix3d inertia = inertiaTensor();
  // Diagonalize inertia tensor V^t * I * V
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(inertia);
  if (eigenSolver.info() != Eigen::Success)
    abort();
  // Rotate to Principal Axes Frame
  this->rotate(eigenSolver.eigenvectors());
  std::cout << eigenSolver.eigenvalues() << std::endl;
}

Molecule & Molecule::operator=(const Molecule & other) {
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
  pointGroup_ = other.pointGroup_;

  return *this;
}

std::ostream & operator<<(std::ostream & os, const Molecule & m) {
  if (m.nAtoms_ != 0) {
    Eigen::IOFormat CleanFmt(6, Eigen::DontAlignCols, " ", "\n", "", "");
    os << "                 Geometry (in Angstrom)" << std::endl;
    os << "   Center            X             Y             Z     \n";
    os << "------------   ------------  ------------  ------------\n";
    for (size_t i = 0; i < m.nAtoms_; ++i) {
      os << std::setw(12) << m.atoms_[i].symbol;
      os << (m.geometry_.col(i).transpose() * bohrToAngstrom()).format(CleanFmt);
      os << std::endl;
    }
    os << "Rotor type: " << rotorTypeList[m.rotor_];
  } else {
    os << "  No atoms in this molecule!" << std::endl;
  }

  return os;
}

Eigen::VectorXd computeMEP(const Molecule & mol, const std::vector<Element> & el) {
  Eigen::VectorXd mep = Eigen::VectorXd::Zero(el.size());
  for (size_t i = 0; i < mol.nAtoms(); ++i) {
    for (size_t j = 0; j < el.size(); ++j) {
      double dist = (mol.geometry().col(i) - el[j].center()).norm();
      mep(j) += mol.charges(i) / dist;
    }
  }
  return mep;
}

Eigen::VectorXd computeMEP(const Molecule & mol, const Eigen::Matrix3Xd & grid) {
  Eigen::VectorXd mep = Eigen::VectorXd::Zero(grid.cols());
  for (size_t i = 0; i < mol.nAtoms(); ++i) {
    for (int j = 0; j < grid.cols(); ++j) {
      double dist = (mol.geometry().col(i) - grid.col(j)).norm();
      mep(j) += mol.charges(i) / dist;
    }
  }
  return mep;
}

Eigen::VectorXd computeMEP(const std::vector<Element> & el,
                           double charge,
                           const Eigen::Vector3d & origin) {
  Eigen::VectorXd mep = Eigen::VectorXd::Zero(el.size());
  for (size_t i = 0; i < el.size(); ++i) {
    double dist = (origin - el[i].center()).norm();
    mep(i) += charge / dist;
  }
  return mep;
}

double GaussEstimate(const Eigen::VectorXd & charges,
                     double permittivity,
                     double correction) {
  return (-charges.sum() * (permittivity - 1) / (permittivity + correction));
}
} // namespace pcm
