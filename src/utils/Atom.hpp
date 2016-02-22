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

#ifndef ATOM_HPP
#define ATOM_HPP

#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

/*! \file Atom.hpp
 *  \struct Atom
 *  \brief A POD describing an atom.
 *  \author Roberto Di Remigio
 *  \date 2011, 2016
 */

struct Atom
{
  /*! Atomic charge */
  double charge;
  /*! Atomic mass */
  double mass;
  /*! Atomic radius */
  double radius;
  /*! Scaling of the atomic radius */
  double radiusScaling;
  /*! Position of the atom */
  Eigen::Vector3d position;
  /*! Name of the element */
  std::string element;
  /*! Atomic symbol */
  std::string symbol;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
  Atom() {}
  Atom(const std::string & elem, const std::string & sym,
      double c, double m, double r,
      const Eigen::Vector3d & coord, double scal = 1.2)
    : charge(c), mass(m), radius(r), radiusScaling(scal),
      position(coord), element(elem), symbol(sym) {}
};

/*! An atom is invalid if it has zero radius */
bool invalid(const Atom & atom);

/*! \brief Returns a reference to a vector<Atom> containing Bondi van der Waals radii.
 *
 * The van der Waals radii are taken from:
 * --- A. Bondi, J. Phys. Chem. 68, 441-451 (1964)
 * complemented with the ones reported in:
 * --- M. Mantina, A. C. Chamberlin, R. Valero, C. J. Cramer, D. G. Truhlar,
 *     J. Phys. Chem. A, 113, 5806-5812 (2009)
 * We are here using Angstrom as in the papers.
 */
std::vector<Atom> & initBondi();

/*! \brief Returns a reference to a vector<Atom> containing UFF radii.
 *
 * The UFF set of radii is taken from:
 * --- A. Rappé, C. J. Casewit, K. S. Colwell, W. A. Goddard, W. M. Skiff,
 *     J. Am. Chem. Soc., 114, 10024-10035 (1992)
 * We are here using Angstrom as in the paper.
 */
std::vector<Atom> & initUFF();

#endif // ATOM_HPP
