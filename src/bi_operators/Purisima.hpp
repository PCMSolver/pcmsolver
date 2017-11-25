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

#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

namespace pcm {
struct BIOperatorData;
namespace cavity {
class Element;
} // namespace cavity
class IGreensFunction;
} // namespace pcm

#include "IBoundaryIntegralOperator.hpp"

namespace pcm {
namespace bi_operators {
/*! \file Purisima.hpp
 *  \class Purisima
 *  \brief Implementation of the double layer operator matrix representation using
 *  one-point collocation and Purisima's strategy for the diagonal of D
 *  \author Roberto Di Remigio
 *  \date 2015, 2016
 *
 *  Calculates the diagonal elements of D as:
 *  \f[
 *  	D_{ii} = -\left(2\pi + \sum_{j\neq i}D_{ij}a_j \right)\frac{1}{a_i}
 *  \f]
 *  The original reference is \cite Purisima1995
 */
class Purisima __final : public IBoundaryIntegralOperator {
public:
  Purisima();
  Purisima(double fac);
  virtual ~Purisima() {}

private:
  /*! Scaling factor for the diagonal elements of the matrix representation of
   * the S operator
   */
  double factor_;
  virtual Eigen::MatrixXd computeS_impl(const std::vector<cavity::Element> & elems,
                                        const IGreensFunction & gf) const __override;
  /*! Computes the matrix representation of the double layer operator by collocation
   *  using the Purisima sum rule to compute the diagonal elements.
   *  \param[in] cav discretized cavity
   *  \param[in] gf  a Green's function
   *
   *  The sum rule for the diagonal elements is:
   *  \f[
   *    D_{ii} = -\left(2\pi + \sum_{j\neq i}D_{ij}a_j \right)\frac{1}{a_i}
   *  \f]
   */
  virtual Eigen::MatrixXd computeD_impl(const std::vector<cavity::Element> & elems,
                                        const IGreensFunction & gf) const __override;
};

IBoundaryIntegralOperator * createPurisima(const BIOperatorData & data);
} // namespace bi_operators
} // namespace pcm
