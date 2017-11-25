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

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "utils/Stencils.hpp"

namespace green {
/*! Returns value of the directional derivative of the function passed for the pair
 * of points p1, p2:
 *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
 * \mathbf{n}_{\mathbf{p}_2}\f$
 *  Notice that this method returns the directional derivative with respect
 *  to the probe point.
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \param[in] functor   function-object handle to the evaluation of the
 * differentiable function
 *  \param[in] normal_p2 the normal vector to p2
 *  \param[in]        p1 first point
 *  \param[in]        p2 second point
 */
template <typename DerivativeTraits>
double derivativeProbe(
    const pcm::function<DerivativeTraits(DerivativeTraits *, DerivativeTraits *)> &
        functor,
    const Eigen::Vector3d & normal_p2,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) {
  DerivativeTraits t1[3], t2[3], der;
  t1[0] = p1(0);
  t1[1] = p1(1);
  t1[2] = p1(2);
  t2[0] = p2(0);
  t2[1] = p2(1);
  t2[2] = p2(2);
  t2[0][1] = normal_p2(0);
  t2[1][1] = normal_p2(1);
  t2[2][1] = normal_p2(2);
  der = functor(t1, t2);
  return der[1];
}

/*! Returns value of the directional derivative of the function passed for the pair
 * of points p1, p2:
 *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
 * \mathbf{n}_{\mathbf{p}_1}\f$
 *  Notice that this method returns the directional derivative with respect
 *  to the source point.
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \param[in] functor   function-object handle to the evaluation of the
 * differentiable function
 *  \param[in] normal_p1 the normal vector to p1
 *  \param[in]        p1 first point
 *  \param[in]        p2 second point
 */
template <typename DerivativeTraits>
double derivativeSource(
    const pcm::function<DerivativeTraits(DerivativeTraits *, DerivativeTraits *)> &
        functor,
    const Eigen::Vector3d & normal_p1,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) {
  DerivativeTraits t1[3], t2[3], der;
  t1[0] = p1(0);
  t1[1] = p1(1);
  t1[2] = p1(2);
  t1[0][1] = normal_p1(0);
  t1[1][1] = normal_p1(1);
  t1[2][1] = normal_p1(2);
  t2[0] = p2(0);
  t2[1] = p2(1);
  t2[2] = p2(2);
  der = functor(t1, t2);
  return der[1];
}

/*! Returns value of the directional derivative of the function passed for the pair
 * of points p1, p2:
 *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
 * \mathbf{n}_{\mathbf{p}_2}\f$
 *  Notice that this method returns the directional derivative with respect
 *  to the probe point.
 *  \param[in] functor   function-object handle to the evaluation of the
 * differentiable function
 *  \param[in] normal_p2 the normal vector to p2
 *  \param[in]        p1 first point
 *  \param[in]        p2 second point
 */
double derivativeProbe(const DifferentiableFunction & functor,
                       const Eigen::Vector3d & normal_p2,
                       const Eigen::Vector3d & p1,
                       const Eigen::Vector3d & p2) {
  return threePointStencil(pcm::bind(functor, pcm::_1, _2), p2, p1, normal_p2);
}

/*! Returns value of the directional derivative of the function passed for the pair
 * of points p1, p2:
 *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
 * \mathbf{n}_{\mathbf{p}_1}\f$
 *  Notice that this method returns the directional derivative with respect
 *  to the source point.
 *  \param[in] functor   function-object handle to the evaluation of the
 * differentiable function
 *  \param[in] normal_p1 the normal vector to p1
 *  \param[in]        p1 first point
 *  \param[in]        p2 second point
 */
double derivativeSource(const DifferentiableFunction & functor,
                        const Eigen::Vector3d & normal_p1,
                        const Eigen::Vector3d & p1,
                        const Eigen::Vector3d & p2) {
  return threePointStencil(pcm::bind(functor, _1, _2), p1, p2, normal_p1);
}

/*! Returns full gradient of the function passed for the pair of points p1, p2:
 *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
 *  Notice that this method returns the gradient with respect to the source point.
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \param[in] functor   function-object handle to the evaluation of the
 * differentiable function
 *  \param[in] p1 first point
 *  \param[in] p2 second point
 */
template <typename DerivativeTraits>
Eigen::Vector3d gradientSource(
    const pcm::function<DerivativeTraits(DerivativeTraits *, DerivativeTraits *)> &
        functor,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) {
  return (Eigen::Vector3d() << derivativeSource(
              functor, Eigen::Vector3d::UnitX(), p1, p2),
          derivativeSource(functor, Eigen::Vector3d::UnitY(), p1, p2),
          derivativeSource(functor, Eigen::Vector3d::UnitZ(), p1, p2))
      .finished();
}

/*! Returns full gradient of the function passed for the pair of points p1, p2:
 *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
 *  Notice that this method returns the gradient with respect to the probe point.
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \param[in] functor   function-object handle to the evaluation of the
 * differentiable function
 *  \param[in] p1 first point
 *  \param[in] p2 second point
 */
template <typename DerivativeTraits>
Eigen::Vector3d gradientProbe(
    const pcm::function<DerivativeTraits(DerivativeTraits *, DerivativeTraits *)> &
        functor,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) {
  return (Eigen::Vector3d() << derivativeProbe(
              functor, Eigen::Vector3d::UnitX(), p1, p2),
          derivativeProbe(functor, Eigen::Vector3d::UnitY(), p1, p2),
          derivativeProbe(functor, Eigen::Vector3d::UnitZ(), p1, p2))
      .finished();
}

/*! Returns full gradient of the function passed for the pair of points p1, p2:
 *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
 *  Notice that this method returns the gradient with respect to the source point.
 *  \param[in] functor   function-object handle to the evaluation of the
 * differentiable function
 *  \param[in] p1 first point
 *  \param[in] p2 second point
 */
Eigen::Vector3d gradientSource(const DifferentiableFunction & functor,
                               const Eigen::Vector3d & p1,
                               const Eigen::Vector3d & p2) {
  return (Eigen::Vector3d() << derivativeSource(
              functor, Eigen::Vector3d::UnitX(), p1, p2),
          derivativeSource(functor, Eigen::Vector3d::UnitY(), p1, p2),
          derivativeSource(functor, Eigen::Vector3d::UnitZ(), p1, p2))
      .finished();
}

/*! Returns full gradient of the function passed for the pair of points p1, p2:
 *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
 *  Notice that this method returns the gradient with respect to the probe point.
 *  \param[in] functor   function-object handle to the evaluation of the
 * differentiable function
 *  \param[in] p1 first point
 *  \param[in] p2 second point
 */
Eigen::Vector3d gradientProbe(const DifferentiableFunction & functor,
                              const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) {
  return (Eigen::Vector3d() << derivativeProbe(
              functor, Eigen::Vector3d::UnitX(), p1, p2),
          derivativeProbe(functor, Eigen::Vector3d::UnitY(), p1, p2),
          derivativeProbe(functor, Eigen::Vector3d::UnitZ(), p1, p2))
      .finished();
}
} // namespace green
