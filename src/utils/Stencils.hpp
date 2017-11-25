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

#include <functional>

#include "Config.hpp"

#include <Eigen/Core>

/*! \typedef DifferentiableFunction
 *  \brief sort of a function pointer to a function of a pair of vectors that can be
 * numerically differentiated
 */
typedef pcm::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)>
    DifferentiableFunction;

/*! \brief Calculate directional derivative using a three-point stencil
 *  \param[in] func function to be differentiated
 *  \param[in] arg_1 first point, the directional derivative is calculated with
 * respect to this point
 *  \param[in] arg_2 second point
 *  \param[in] direction the direction in which the directional derivative is to be
 * evaluated
 *  \param[in] step finite difference value for the stencil
 */
inline double threePointStencil(const DifferentiableFunction & func,
                                const Eigen::Vector3d & arg_1,
                                const Eigen::Vector3d & arg_2,
                                const Eigen::Vector3d & direction,
                                double step = 1.0e-04) {
  // f(x-h)
  Eigen::Vector3d delta_m1 = arg_1 - direction * step / direction.norm();
  // f(x+h)
  Eigen::Vector3d delta_1 = arg_1 + direction * step / direction.norm();

  Eigen::Vector2d stencil;
  stencil << -0.5, 0.5;
  Eigen::Vector2d function_values;
  function_values << func(delta_m1, arg_2), func(delta_1, arg_2);

  return (function_values.dot(stencil) / step);
}

/*! \brief Calculate directional derivative using a five-point stencil
 *  \param[in] func function to be differentiated
 *  \param[in] arg_1 first point, the directional derivative is calculated with
 * respect to this point
 *  \param[in] arg_2 second point
 *  \param[in] direction the direction in which the directional derivative is to be
 * evaluated
 *  \param[in] step finite difference value for the stencil
 */
inline double fivePointStencil(const DifferentiableFunction & func,
                               const Eigen::Vector3d & arg_1,
                               const Eigen::Vector3d & arg_2,
                               const Eigen::Vector3d & direction,
                               double step = 1.0e-04) {
  // f(x-2h)
  Eigen::Vector3d delta_m2 = arg_1 - 2.0 * direction * step / direction.norm();
  // f(x-h)
  Eigen::Vector3d delta_m1 = arg_1 - direction * step / direction.norm();
  // f(x+h)
  Eigen::Vector3d delta_1 = arg_1 + direction * step / direction.norm();
  // f(x+2h)
  Eigen::Vector3d delta_2 = arg_1 + 2.0 * direction * step / direction.norm();

  Eigen::Vector4d stencil;
  stencil << 1.0 / 12.0, -2.0 / 3.0, 2.0 / 3.0, -1.0 / 12.0;
  Eigen::Vector4d function_values;
  function_values << func(delta_m2, arg_2), func(delta_m1, arg_2),
      func(delta_1, arg_2), func(delta_2, arg_2);

  return (function_values.dot(stencil) / step);
}

/*! \brief Calculate directional derivative using a seven-point stencil
 *  \param[in] func function to be differentiated
 *  \param[in] arg_1 first point, the directional derivative is calculated with
 * respect to this point
 *  \param[in] arg_2 second point
 *  \param[in] direction the direction in which the directional derivative is to be
 * evaluated
 *  \param[in] step finite difference value for the stencil
 */
inline double sevenPointStencil(const DifferentiableFunction & func,
                                const Eigen::Vector3d & arg_1,
                                const Eigen::Vector3d & arg_2,
                                const Eigen::Vector3d & direction,
                                double step = 1.0e-04) {
  // f(x-3h)
  Eigen::Vector3d delta_m3 = arg_1 - 3.0 * direction * step / direction.norm();
  // f(x-2h)
  Eigen::Vector3d delta_m2 = arg_1 - 2.0 * direction * step / direction.norm();
  // f(x-h)
  Eigen::Vector3d delta_m1 = arg_1 - direction * step / direction.norm();
  // f(x+h)
  Eigen::Vector3d delta_1 = arg_1 + direction * step / direction.norm();
  // f(x+2h)
  Eigen::Vector3d delta_2 = arg_1 + 2.0 * direction * step / direction.norm();
  // f(x+3h)
  Eigen::Vector3d delta_3 = arg_1 + 3.0 * direction * step / direction.norm();

  Eigen::Matrix<double, 6, 1> stencil;
  stencil << -1.0 / 60.0, 3.0 / 20.0, -3.0 / 4.0, 3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0;
  Eigen::Matrix<double, 6, 1> function_values;
  function_values << func(delta_m3, arg_2), func(delta_m2, arg_2),
      func(delta_m1, arg_2), func(delta_1, arg_2), func(delta_2, arg_2),
      func(delta_3, arg_2);

  return (function_values.dot(stencil) / step);
}
