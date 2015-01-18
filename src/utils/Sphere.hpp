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

#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

/*! \file Sphere.hpp
 *  \class Sphere
 *  \brief Class describing a sphere.
 *  \author Roberto Di Remigio
 *  \date 2011
 */

class Sphere
{
public:
    Sphere() {}
    Sphere(const Eigen::Vector3d & center, double radius,
           const std::string & colour = "Violet" )
        : center_(center), radius_(radius), colour_(colour) {}
    ~Sphere() {}
    double radius() const { return radius_; }
    void radius(double r) { radius_ = r; }
    const Eigen::Vector3d & center() const { return center_; }
    double center(int i) const { return center_(i); }
    void center(Eigen::Vector3d & coord) { center_ = coord; }
    const std::string & colour() const { return colour_; }
    void colour(std::string & col) { colour_ = col; }
    friend inline void swap(Sphere & left, Sphere & right);
    inline void swap(Sphere & other);
    /// Assignment operator.
    Sphere & operator=(Sphere other);
    friend std::ostream& operator<<(std::ostream & os, Sphere & sph) {
        return sph.printObject(os);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    Eigen::Vector3d center_;
    double radius_;
    std::string colour_;
    std::ostream & printObject(std::ostream & os);
};

/*! \fn inline void transfer_spheres(const std::vector<Sphere> & spheres, Eigen::Matrix3Xd & sphereCenter, Eigen::VectorXd & sphereRadius)
 *  \brief Transfer info from std::vector<Sphere> to Eigen objects.
 *  \param[in] spheres list of spheres as std::vector<Sphere>
 *  \param[out] sphereCenter sphere centers as Eigen::Matrix3Xd (xyz * nSpheres)
 *  \param[out] sphereRadius sphere radii as Eigen::VectorXd
 *
 *  This is used in the Cavity.hpp constructor
 */
inline void transfer_spheres(const std::vector<Sphere> & spheres, Eigen::Matrix3Xd & sphereCenter, Eigen::VectorXd & sphereRadius)
{
    size_t nSpheres = spheres.size();
    sphereCenter.resize(Eigen::NoChange, nSpheres);
    sphereRadius.resize(nSpheres);
    for (size_t i = 0; i < nSpheres; ++i) {
        sphereCenter.col(i) = spheres[i].center();
        sphereRadius(i) = spheres[i].radius();
    }
}

#endif // SPHERE_HPP
