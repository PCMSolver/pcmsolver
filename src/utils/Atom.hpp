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

#ifndef ATOM_HPP
#define ATOM_HPP

#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

/*! \file Atom.hpp
 *  \class Atom
 *  \brief Class describing an atom.
 *  \author Roberto Di Remigio
 *  \date 2011
 *
 *  This class contains all the radii sets available in the module.
 *  They can be obtained by the client through static functions initializing
 *  the sets via the lazy evaluation idiom.
 */

class Atom
{
public:
    Atom() {}
    Atom(const std::string & element, const std::string & symbol, double charge,
         double radius, const Eigen::Vector3d & coord, double scaling = 1.0,
         const std::string & colour = "Violet")
        : atomElement_(element), atomSymbol_(symbol), atomCharge_(charge),
          atomRadius_(radius),
          atomCoord_(coord), atomRadiusScaling_(scaling), atomColour_(colour) {}
    Atom(const std::string & element, const std::string & symbol, double charge,
         double radius)
        : atomElement_(element), atomSymbol_(symbol), atomCharge_(charge),
          atomRadius_(radius) {
        Eigen::Vector3d Origin(0.0, 0.0, 0.0);
        std::string colour = "Violet";
        atomCoord_ = Origin;
        atomColour_ = colour;
        atomRadiusScaling_ = 1.0;
    }
    ~Atom() {}
    std::string atomElement() const { return atomElement_; }
    void atomElement(const std::string & element) { atomElement_ = element; }
    std::string atomSymbol() const { return atomSymbol_; }
    void atomSymbol(const std::string & symbol) { atomSymbol_ = symbol; }
    const Eigen::Vector3d & atomCoord() const { return atomCoord_; }
    void atomCoord(Eigen::Vector3d & coord) { atomCoord_ = coord; }
    double atomCharge() const { return atomCharge_; }
    void atomCharge(double charge) { atomCharge_ = charge; }
    /*! \brief Returns the atomic radius in Angstrom
     */
    double atomRadius() const { return atomRadius_; }
    void atomRadius(double radius) { atomRadius_ = radius; }
    double atomRadiusScaling() const { return atomRadiusScaling_; }
    void atomRadiusScaling(double scaling) { atomRadiusScaling_ = scaling; }
    std::string atomColour() const { return atomColour_; }
    void atomColour(const std::string & colour) { atomColour_ = colour; }
    /*! \brief Returns a reference to a vector<Atom> containing Bondi van der Waals radii.
     *
     * The van der Waals radii are taken from:
     * --- A. Bondi, J. Phys. Chem. 68, 441-451 (1964)
     * complemented with the ones reported in:
     * --- M. Mantina, A. C. Chamberlin, R. Valero, C. J. Cramer, D. G. Truhlar,
     *     J. Phys. Chem. A, 113, 5806-5812 (2009)
     * We are here using Angstrom as in the papers.
     */
    static std::vector<Atom> & initBondi();
    /*! \brief Returns a reference to a vector<Atom> containing UFF radii.
     *
     * The UFF set of radii is taken from:
     * --- A. Rappé, C. J. Casewit, K. S. Colwell, W. A. Goddard, W. M. Skiff,
     *     J. Am. Chem. Soc., 114, 10024-10035 (1992)
     * We are here using Angstrom as in the paper.
     */
    static std::vector<Atom> & initUFF();
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    std::string atomElement_;
    std::string atomSymbol_;
    double atomCharge_;
    double atomRadius_;
    Eigen::Vector3d atomCoord_;
    double atomRadiusScaling_;
    std::string atomColour_;
};

#endif // ATOM_HPP
