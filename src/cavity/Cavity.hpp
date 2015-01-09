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

#ifndef CAVITY_HPP
#define CAVITY_HPP

#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Molecule.hpp"
#include "Sphere.hpp"
#include "Symmetry.hpp"

/*!
 * \file Cavity.hpp
 * \class Cavity
 * \brief Abstract Base Class for cavities.
 * \author Krzysztof Mozgawa
 * \date 2011
 * 
 * This class represents a cavity made of spheres, its surface being discretized in
 * terms of finite elements.
 */

class Cavity
{
protected:
    std::vector<Sphere> spheres_;
    Molecule molecule_;
    int nElements_;
    int nIrrElements_;
    bool built;
    Eigen::Matrix3Xd elementCenter_;
    Eigen::Matrix3Xd elementNormal_;
    Eigen::VectorXd elementArea_;
    int nSpheres_;
    Eigen::Matrix3Xd elementSphereCenter_;
    Eigen::VectorXd elementRadius_;
    Eigen::Matrix3Xd sphereCenter_;
    Eigen::VectorXd sphereRadius_;
private:
    /*! \brief Creates the cavity and discretizes its surface.
     *
     *  Has to be implemented by classes lower down in the inheritance hierarchy
     */
    virtual void makeCavity() = 0;
    virtual std::ostream & printCavity(std::ostream & os) = 0;
public:
    //! Default constructor
    Cavity() : nElements_(0), built(false) {}
    /*! \brief Constructor from list of spheres 
     *  \param[in] sph the list of spheres
     *
     *  Only used when we have to deal with a single sphere, i.e. in the unit tests
     */
    Cavity(const std::vector<Sphere> & sph) : spheres_(sph), built(false) {
	molecule_ = Molecule(spheres_); 
        nSpheres_ = spheres_.size();
	transfer_spheres(spheres_, sphereCenter_, sphereRadius_);
    }
    /*! \brief Constructor from Molecule
     *  \param[in] molec the molecular aggregate
     */
    Cavity(const Molecule & molec) : spheres_(molec.spheres()), molecule_(molec), built(false) {
        nSpheres_ = spheres_.size();
	transfer_spheres(spheres_, sphereCenter_, sphereRadius_);
    }
    virtual ~Cavity() {}
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
    int size() { return nElements_; }
    int size() const { return nElements_; }
    int irreducible_size() { return nIrrElements_; }
    int irreducible_size() const { return nIrrElements_; }
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
    const Eigen::Matrix3Xd & elementSphereCenter() const { return elementSphereCenter_; }
    bool isBuilt() { return built; }
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
    friend std::ostream & operator<<(std::ostream & os, Cavity & cavity) {
        return cavity.printCavity(os);
    }
};

#endif // CAVITY_HPP
