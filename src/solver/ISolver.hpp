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

#include "Config.hpp"

/*! \file ISolver.hpp */

class ICavity;
class IGreensFunction;

namespace pcm {
class IBoundaryIntegralOperator;
} // namespace pcm

namespace pcm {
/*! \class ISolver
 *  \brief Abstract Base Class for solvers inheritance hierarchy.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2011, 2015, 2016
 *
 *  We use the Non-Virtual Interface idiom.
 */
class ISolver {
public:
  ISolver() : built_(false), isotropic_(true) {}
  virtual ~ISolver() {}

  /*! \brief Calculation of the PCM matrix
   *  \param[in] cavity the cavity to be used
   *  \param[in] gf_i Green's function inside the cavity
   *  \param[in] gf_o Green's function outside the cavity
   *  \param[in] op integrator strategy for the single and double layer operators
   */
  void buildSystemMatrix(const ICavity & cavity,
                         const IGreensFunction & gf_i,
                         const IGreensFunction & gf_o,
                         const IBoundaryIntegralOperator & op) {
    buildSystemMatrix_impl(cavity, gf_i, gf_o, op);
  }
  /*! \brief Returns the ASC given the MEP and the desired irreducible representation
   *  \param[in] potential the vector containing the MEP at cavity points
   *  \param[in] irrep the irreducible representation of the MEP and ASC
   */
  Eigen::VectorXd computeCharge(const Eigen::VectorXd & potential,
                                int irrep = 0) const {
    if (!built_)
      PCMSOLVER_ERROR("PCM matrix not calculated yet");
    return computeCharge_impl(potential, irrep);
  }

  friend std::ostream & operator<<(std::ostream & os, ISolver & solver) {
    return solver.printSolver(os);
  }

protected:
  /*! Whether the system matrix has been built */
  bool built_;
  /*! Whether the solver is isotropic */
  bool isotropic_;

  /*! \brief Calculation of the PCM matrix
   *  \param[in] cavity the cavity to be used
   *  \param[in] gf_i Green's function inside the cavity
   *  \param[in] gf_o Green's function outside the cavity
   *  \param[in] op integrator strategy for the single and double layer operators
   */
  virtual void buildSystemMatrix_impl(const ICavity & cavity,
                                      const IGreensFunction & gf_i,
                                      const IGreensFunction & gf_o,
                                      const IBoundaryIntegralOperator & op) = 0;
  /*! \brief Returns the ASC given the MEP and the desired irreducible representation
   *  \param[in] potential the vector containing the MEP at cavity points
   *  \param[in] irrep the irreducible representation of the MEP and ASC
   */
  virtual Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential,
                                             int irrep = 0) const = 0;
  virtual std::ostream & printSolver(std::ostream & os) = 0;
};
} // namespace pcm
