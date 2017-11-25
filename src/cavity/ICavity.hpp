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

#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "Element.hpp"
#include "utils/Molecule.hpp"
#include "utils/Sphere.hpp"
#include "utils/Symmetry.hpp"

/*! \file ICavity.hpp */

namespace pcm {
using cavity::Element;
using utils::Sphere;

/*! \class ICavity
 * \brief Abstract Base Class for cavities.
 * \author Krzysztof Mozgawa
 * \date 2011
 *
 * This class represents a cavity made of spheres, its surface being discretized in
 * terms of finite elements.
 */
class ICavity {
protected:
  /// List of spheres
  std::vector<Sphere> spheres_;
  /// The molecule to be wrapped by the cavity
  Molecule molecule_;
  /// Number of finite elements generated
  PCMSolverIndex nElements_;
  /// Number of irreducible finite elements
  PCMSolverIndex nIrrElements_;
  /// Whether the cavity has been built
  bool built;
  /// Coordinates of elements centers
  Eigen::Matrix3Xd elementCenter_;
  /// Outward-pointing normal vectors to the elements centers
  Eigen::Matrix3Xd elementNormal_;
  /// Elements areas
  Eigen::VectorXd elementArea_;
  /// Number of spheres
  int nSpheres_;
  /// Centers of the sphere the elements belong to
  Eigen::Matrix3Xd elementSphereCenter_;
  /// Radii of the sphere the elements belong to
  Eigen::VectorXd elementRadius_;
  /// Spheres centers
  Eigen::Matrix3Xd sphereCenter_;
  /// Spheres radii
  Eigen::VectorXd sphereRadius_;
  /// List of finite elements
  std::vector<Element> elements_;
  /// Molecular point group
  Symmetry pointGroup_;

private:
  /*! \brief Creates the cavity and discretizes its surface.
   *
   *  Has to be implemented by classes lower down in the inheritance hierarchy
   */
  virtual void makeCavity() = 0;
  virtual std::ostream & printCavity(std::ostream & os) = 0;

public:
  //! Default constructor
  ICavity();
  /*! \brief Constructor from a single sphere
   *  \param[in] sph the sphere
   *
   *  Only used when we have to deal with a single sphere, i.e. in the unit tests
   */
  ICavity(const Sphere & sph);
  /*! \brief Constructor from list of spheres
   *  \param[in] sph the list of spheres
   *
   *  Only used when we have to deal with a single sphere, i.e. in the unit tests
   */
  ICavity(const std::vector<Sphere> & sph);
  /*! \brief Constructor from Molecule
   *  \param[in] molec the molecular aggregate
   */
  ICavity(const Molecule & molec);
  virtual ~ICavity() {}
  Eigen::Matrix3Xd & elementCenter() { return elementCenter_; }
  const Eigen::Matrix3Xd & elementCenter() const { return elementCenter_; }
  Eigen::Vector3d elementCenter(int i) { return elementCenter_.col(i); }
  Eigen::Vector3d elementCenter(int i) const { return elementCenter_.col(i); }
  Eigen::Matrix3Xd & elementNormal() { return elementNormal_; }
  const Eigen::Matrix3Xd & elementNormal() const { return elementNormal_; }
  Eigen::Vector3d elementNormal(int i) { return elementNormal_.col(i); }
  Eigen::Vector3d elementNormal(int i) const { return elementNormal_.col(i); }
  Eigen::VectorXd & elementArea() { return elementArea_; }
  const Eigen::VectorXd & elementArea() const { return elementArea_; }
  double elementArea(int i) { return elementArea_(i); }
  double elementArea(int i) const { return elementArea_(i); }
  PCMSolverIndex size() { return nElements_; }
  PCMSolverIndex size() const { return nElements_; }
  PCMSolverIndex irreducible_size() { return nIrrElements_; }
  PCMSolverIndex irreducible_size() const { return nIrrElements_; }
  Symmetry pointGroup() const { return molecule_.pointGroup(); }
  std::vector<Sphere> & spheres() { return spheres_; }
  const std::vector<Sphere> & spheres() const { return spheres_; }
  int nSpheres() { return nSpheres_; }
  int nSpheres() const { return nSpheres_; }
  Eigen::VectorXd & sphereRadius() { return sphereRadius_; }
  const Eigen::VectorXd & sphereRadius() const { return sphereRadius_; }
  Eigen::Matrix3Xd & sphereCenter() { return sphereCenter_; }
  const Eigen::Matrix3Xd & sphereCenter() const { return sphereCenter_; }
  Eigen::VectorXd & elementRadius() { return elementRadius_; }
  const Eigen::VectorXd & elementRadius() const { return elementRadius_; }
  double elementRadius(int i) { return elementRadius_(i); }
  double elementRadius(int i) const { return elementRadius_(i); }
  Eigen::Matrix3Xd & elementSphereCenter() { return elementSphereCenter_; }
  const Eigen::Matrix3Xd & elementSphereCenter() const {
    return elementSphereCenter_;
  }
  bool isBuilt() { return built; }
  const std::vector<Element> & elements() const { return elements_; }
  const Element & elements(int i) const { return elements_[i]; }
  /*! \brief Save cavity specification to file.
   *
   *  The cavity specification contains:
   *   0. the number of finite elements, nElements;
   *   1. the weight of the finite elements, elementArea;
   *   2. the radius of the finite elements, elementRadius;
   *   3. the centers of the finite elements, elementCenter;
   *   4. the normal vectors relative to the centers, elementNormal.
   *  Each of these objects is saved in a separate .npy binary file
   *  and compressed into one .npz file.
   *  Notice that this is just the minimal set of data needed to
   *  restart an energy calculation.
   */
  virtual void saveCavity(const std::string & fname = "cavity.npz");
  /*! \brief Load cavity specification from file.
  */
  virtual void loadCavity(const std::string & fname = "cavity.npz");
  friend std::ostream & operator<<(std::ostream & os, ICavity & cavity) {
    return cavity.printCavity(os);
  }
};
} // namespace pcm
