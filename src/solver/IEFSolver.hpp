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
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

/*! \file IEFSolver.hpp */

namespace pcm {
class ICavity;
class IGreensFunction;
class IBoundaryIntegralOperator;
struct SolverData;
} // namespace pcm

#include "ISolver.hpp"

namespace pcm {
namespace solver {
/*! \class IEFSolver
 *  \brief IEFPCM, collocation-based solver
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2011, 2015, 2016
 *
 *  \note We store the non-Hermitian, symmetry-blocked T(epsilon) and Rinfinity
 *matrices.
 *  The ASC is obtained by multiplying the MEP by Rinfinity and then using a
 *partially
 *  pivoted LU decomposition of T(epsilon) on the resulting vector.
 *  In case the polarization weights are requested, we use the approach suggested in
 *  \cite Barone2004-ae.
 *  First, the adjoint problem is solved:
 *  \f[
 *      \mathbf{T}_\varepsilon^\dagger \tilde{v} = v
 *  \f]
 *  Also in this case a partially pivoted LU decomposition is used.
 *  The "transposed" ASC is obtained by the matrix-vector multiplication:
 *  \f[
 *      q^* = \mathbf{R}_\infty^\dagger \tilde{v}
 *  \f]
 *  Eventually, the two sets of charges are summed and divided by 2
 *  This avoids computing and storing the inverse explicitly, at the expense of
 *storing
 *  both T(epsilon) and Rinfinity.
 */
class IEFSolver : public ISolver {
public:
  IEFSolver() {}
  /*! \brief Construct solver
   *  \param[in] symm whether the system matrix has to be symmetrized
   */
  IEFSolver(bool symm) : ISolver(), hermitivitize_(symm) {}
  virtual ~IEFSolver() {}
  /*! \brief Builds PCM matrix for an anisotropic environment
   *  \param[in] cavity the cavity to be used.
   *  \param[in] gf_i Green's function inside the cavity
   *  \param[in] gf_o Green's function outside the cavity
   *  \param[in] op integrator strategy for the single and double layer operators
   */
  void buildAnisotropicMatrix(const ICavity & cavity,
                              const IGreensFunction & gf_i,
                              const IGreensFunction & gf_o,
                              const IBoundaryIntegralOperator & op);
  /*! \brief Builds PCM matrix for an isotropic environment
   *  \param[in] cavity the cavity to be used.
   *  \param[in] gf_i Green's function inside the cavity
   *  \param[in] gf_o Green's function outside the cavity
   *  \param[in] op integrator strategy for the single and double layer operators
   */
  void buildIsotropicMatrix(const ICavity & cavity,
                            const IGreensFunction & gf_i,
                            const IGreensFunction & gf_o,
                            const IBoundaryIntegralOperator & op);
  friend std::ostream & operator<<(std::ostream & os, IEFSolver & solver) {
    return solver.printSolver(os);
  }

private:
  /*! Whether the system matrix has to be symmetrized */
  bool hermitivitize_;
  /*! T(epsilon) matrix, not symmetry blocked */
  Eigen::MatrixXd Tepsilon_;
  /*! T(epsilon) matrix, symmetry blocked form */
  std::vector<Eigen::MatrixXd> blockTepsilon_;
  /*! R_infinity matrix, not symmetry blocked */
  Eigen::MatrixXd Rinfinity_;
  /*! R_infinity matrix, symmetry blocked form */
  std::vector<Eigen::MatrixXd> blockRinfinity_;

  virtual void buildSystemMatrix_impl(const ICavity & cavity,
                                      const IGreensFunction & gf_i,
                                      const IGreensFunction & gf_o,
                                      const IBoundaryIntegralOperator & op)
      __override;
  virtual Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential,
                                             int irrep = 0) const __override;
  virtual std::ostream & printSolver(std::ostream & os) __override;
};

ISolver * createIEFSolver(const SolverData & data);
} // namespace solver
} // namespace pcm
