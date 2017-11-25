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
/*! \file Numerical.hpp
 *  \class Numerical
 *  \brief Implementation of the single and double layer operators matrix
 *representation using one-point collocation
 *  \author Roberto Di Remigio
 *  \date 2015, 2016
 *
 *  Calculates the diagonal elements of S and D by collocation, using numerical
 *integration.
 */
class Numerical __final : public IBoundaryIntegralOperator {
private:
  virtual Eigen::MatrixXd computeS_impl(const std::vector<cavity::Element> & elems,
                                        const IGreensFunction & gf) const __override;
  virtual Eigen::MatrixXd computeD_impl(const std::vector<cavity::Element> & elems,
                                        const IGreensFunction & gf) const __override;
};

IBoundaryIntegralOperator * createNumerical(const BIOperatorData & /* data */);

/*! \typedef KernelS
 *  \brief functor handle to the kernelS method in IGreensFunction
 */
typedef pcm::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)>
    KernelS;

/*! \typedef KernelD
 *  \brief functor handle to the kernelD method in IGreensFunction
 */
typedef pcm::function<double(const Eigen::Vector3d &,
                             const Eigen::Vector3d &,
                             const Eigen::Vector3d &)>
    KernelD;

/*! \brief Integrates a single layer type operator on a single spherical polygon
 *  \date 2014
 *  \tparam PhiPoints Gaussian rule to be used in the angular phi integration
 *  \tparam ThetaPoints Gaussian rule to be used in the angular theta integration
 *
 *  This is needed for the numerical evaluation of the diagonal elements when using
 *  centroid collocation.
 */
template <int PhiPoints, int ThetaPoints>
double integrateS(const KernelS & F, const cavity::Element & e);

/*! \brief Integrates a double layer type operator on a single spherical polygon
 *  \date 2014
 *  \tparam PhiPoints Gaussian rule to be used in the angular phi integration
 *  \tparam ThetaPoints Gaussian rule to be used in the angular theta integration
 *
 *  This is needed for the numerical evaluation of the diagonal elements when using
 *  centroid collocation.
 */
template <int PhiPoints, int ThetaPoints>
double integrateD(const KernelD & F, const cavity::Element & e);
} // namespace bi_operators
} // namespace pcm
