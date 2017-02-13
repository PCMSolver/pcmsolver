/**
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

#ifndef ATOM_HPP
#define ATOM_HPP

#include <map>
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

namespace pcm {
namespace utils {
struct Atom {
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
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See
                                     http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
                                     */
      Atom()
      : charge(0.0),
        mass(0.0),
        radius(0.0),
        radiusScaling(0.0),
        position(Eigen::Vector3d::Zero()),
        element("Dummy"),
        symbol("Du") {}
  Atom(const std::string & elem,
       const std::string & sym,
       double c,
       double m,
       double r,
       const Eigen::Vector3d & coord,
       double scal = 1.0)
      : charge(c),
        mass(m),
        radius(r),
        radiusScaling(scal),
        position(coord),
        element(elem),
        symbol(sym) {}
};

typedef tuple<std::string, std::vector<Atom> > RadiiSet;

namespace detail {
/*! \brief Returns a vector<Atom> containing Bondi van der Waals
 *radii.
 *
 * The van der Waals radii are taken from:
 * --- A. Bondi, J. Phys. Chem. 68, 441-451 (1964)
 * complemented with the ones reported in:
 * --- M. Mantina, A. C. Chamberlin, R. Valero, C. J. Cramer, D. G. Truhlar,
 *     J. Phys. Chem. A, 113, 5806-5812 (2009)
 * We are here using Angstrom as in the papers.
 */
RadiiSet initBondi();

/*! \brief Returns a vector<Atom> containing UFF radii.
 *
 * The UFF set of radii is taken from:
 * --- A. Rappé, C. J. Casewit, K. S. Colwell, W. A. Goddard, W. M. Skiff,
 *     J. Am. Chem. Soc., 114, 10024-10035 (1992)
 * We are here using Angstrom as in the paper.
 */
RadiiSet initUFF();

/*! \brief Returns a vector<Atom> containing Allinger's MM3 radii.
 *
 * The MM3 set of radii is taken from:
 * --- N. L. Allinger, X. Zhou, J. Bergsma,
 *     J. Mol. Struct. (THEOCHEM), 312, 69-83 (1994)
 * We are here using Angstrom as in the paper.
 *
 * \note We *divide* the values reported in the paper by 1.2, as done in
 * the ADF program package.
 */
RadiiSet initAllinger();
} // namespace detail

class Factory __final {
public:
  /*! \brief Callback function for object creation
   *  Returns a std::vector<Atom> by reference
   */
  typedef function<RadiiSet()> Create;

private:
  /*! std::map from the object type identifier (a string) to its callback function */
  typedef std::map<std::string, Create> CallbackMap;
  /*! std::pair of an object type identifier and a callback function */
  typedef typename CallbackMap::value_type CallbackPair;
  /*! const iterator */
  typedef typename CallbackMap::const_iterator CallbackConstIter;

public:
  /*! \brief Returns true on successful registration of the objID
   * \param[in] objID  the object's identification string
   * \param[in] functor the creation function related to the object type given
   */
  bool registerObject(const std::string & objID, const Create & functor) {
    return callbacks_.insert(CallbackPair(objID, functor)).second;
  }
  /*! \brief Returns true if objID was already registered
   *  \param objID the object's identification string
   */
  bool unRegisterObject(const std::string & objID) {
    return callbacks_.erase(objID) == 1;
  }
  /*! \brief Calls the appropriate creation functor, based on the passed objID
   *  \param[in] objID the object's identification string
   *  \param[in] data  input data for the creation of the object
   */
  RadiiSet create(const std::string & objID) {
    if (objID.empty())
      PCMSOLVER_ERROR("No object identification string provided to the Factory.");
    CallbackConstIter i = callbacks_.find(objID);
    if (i == callbacks_.end())
      PCMSOLVER_ERROR("The unknown object ID " + objID +
                      " occurred in the Factory.");
    return (i->second)();
  }

  Factory() {}
  ~Factory() { callbacks_.clear(); }

private:
  CallbackMap callbacks_;
};

inline Factory bootstrapRadiiSet() {
  Factory factory_;

  const bool bondi = factory_.registerObject("BONDI", detail::initBondi);
  if (!bondi)
    PCMSOLVER_ERROR("Subscription of Bondi radii set to factory failed!");

  const bool uff = factory_.registerObject("UFF", detail::initUFF);
  if (!uff)
    PCMSOLVER_ERROR("Subscription of UFF radii set to factory failed!");

  const bool allinger = factory_.registerObject("ALLINGER", detail::initAllinger);
  if (!allinger)
    PCMSOLVER_ERROR("Subscription of Allinger's MM3 radii set to factory failed!");

  return factory_;
}
} // namespace utils

/*! An atom is invalid if it has zero radius */
bool invalid(const utils::Atom & atom);
} // namespace pcm

#endif // ATOM_HPP
