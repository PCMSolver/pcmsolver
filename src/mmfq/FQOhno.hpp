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

#ifndef FQOHNO_HPP
#define FQOHNO_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

#include "utils/MMFQ.hpp"

/*! \file FQOhno.hpp
 *  \class FQOhno
 *  \brief Solver for MMFQ with Ohno interaction kernel
 *  \author Roberto Di Remigio
 *  \date 2017
 *
 *  \note We store the D_\lambda matrix and use a robust Cholesky decomposition
 *  to solve for the charges.
 *  This avoids computing and storing the inverse explicitly.
 */

namespace pcm {
namespace mmfq {
class FQOhno __final {
public:
  FQOhno() {}
  /*! \brief Construct solver
   *  \param[in] the fluctuating charges force field
   */
  FQOhno(const utils::MMFQ & ff);

  ~FQOhno() {}
  friend std::ostream & operator<<(std::ostream & os, FQOhno & solver) {
    return solver.printSolver(os);
  }
  Eigen::VectorXd computeCharge(const Eigen::VectorXd & potential) const {
    return computeCharge_impl(potential);
  }

private:
  bool built_;
  utils::MMFQ mmfq_;
  /*! D_\lambda matrix, not symmetry blocked */
  Eigen::MatrixXd Dlambda_;

  void buildSystemMatrix_impl();
  Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential) const;
  std::ostream & printSolver(std::ostream & os);
};
} // namespace mmfq
} // namespace pcm

#endif // FQOHNO_HPP
