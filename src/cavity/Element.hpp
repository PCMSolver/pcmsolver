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

#include <ostream>
#include <vector>

#include "Config.hpp"

#include "utils/Sphere.hpp"

/*! \file Element.hpp
 *  \class Element
 *  \brief Element data structure
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  Data structure containing relevant information about a finite element
 *  making up the cavity
 */

namespace pcm {
using utils::Sphere;
namespace cavity {
class Element __final {
public:
  Element(int nv,
          int isphe,
          double w,
          const Eigen::Vector3d & c,
          const Eigen::Vector3d & n,
          bool i,
          const Sphere & s,
          const Eigen::Matrix3Xd & v,
          const Eigen::Matrix3Xd & a)
      : nVertices_(nv),
        iSphere_(isphe),
        area_(w),
        center_(c),
        normal_(n),
        irreducible_(i),
        sphere_(s),
        vertices_(v),
        arcs_(a) {}
  ~Element() {}

  int nVertices() const { return nVertices_; }
  int iSphere() const { return iSphere_; }
  double area() const { return area_; }
  Eigen::Vector3d center() const { return center_; }
  Eigen::Vector3d normal() const { return normal_; }
  bool irreducible() const { return irreducible_; }
  Sphere sphere() const { return sphere_; }
  Eigen::Matrix3Xd vertices() const { return vertices_; }
  Eigen::Matrix3Xd arcs() const { return arcs_; }
  /*! \brief Calculate arrays of polar and azimuthal angles for element vertices
   *  \param[in] t_ tangent vector
   *  \param[in] b_ bitangent vector
   *  \param[out] theta   polar angles of vertices
   *  \param[out] phi     azimuthal angles of vertices
   *  \param[out] phinumb contains azimuth for vertices of subtriangles
   *  \param[out] numb    contains index of polar angle for vertices of subtriangles
   */
  void spherical_polygon(Eigen::Vector3d & t_,
                         Eigen::Vector3d & b_,
                         std::vector<double> & theta,
                         std::vector<double> & phi,
                         std::vector<double> & phinumb,
                         std::vector<int> & numb) const;

  friend std::ostream & operator<<(std::ostream & os, Element & element) {
    return element.printElement(os);
  }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
  /// Number of vertices of the finite element
  int nVertices_;
  /// Index of the sphere the finite element belongs to
  int iSphere_;
  /// Area of the finite element
  double area_;
  /// Center of the finite element
  Eigen::Vector3d center_;
  /// Outward-pointing normal to the center of the finite element
  Eigen::Vector3d normal_;
  /// Whether the finite element is irreducible
  bool irreducible_;
  /// The sphere the finite element belongs to
  Sphere sphere_;
  /// Coordinates of the finite element vertices (dimension 3*nVertices_)
  Eigen::Matrix3Xd vertices_;
  /// Coordinates of the centers of the arcs defining the edges of the finite element
  /// (dimension 3*nVertices_)
  Eigen::Matrix3Xd arcs_;
  virtual std::ostream & printElement(std::ostream & os) {
    os << "Finite element" << std::endl;
    os << "Center\n" << center_.transpose() << std::endl;
    os << "Normal\n" << normal_.transpose() << std::endl;
    os << "Weight " << area_ << std::endl;
    os << sphere_ << std::endl;
    os << "Number of vertices and arcs " << nVertices_ << std::endl;
    os << "Vertices\n" << vertices_ << std::endl;
    os << "Arcs\n" << arcs_;
    return os;
  }
};

/*! \brief Calculate tangent and bitangent vector for the representative point
 *  \param[in]  n_ normal vector
 *  \param[out] t_ tangent vector
 *  \param[out] b_ bitangent vector
 */
namespace detail {
void tangent_and_bitangent(const Eigen::Vector3d & n_,
                           Eigen::Vector3d & t_,
                           Eigen::Vector3d & b_);
} // namespace detail
} // namespace cavity
} // namespace pcm
