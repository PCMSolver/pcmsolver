/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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

#include <limits>

#include "Config.hpp"

/*! \file RungeKutta.hpp */

namespace pcm {
namespace utils {
namespace detail {
/*! \return t1 < t2 if dt > 0 and t1 > t2 if dt < 0 with epsilon accuracy
 *  \note This function is part of Boost.Odeint
 */
template <typename T> bool less_with_sign(T t1, T t2, T dt) {
  if (dt > 0) { // return t1 < t2
    return t2 - t1 > std::numeric_limits<T>::epsilon();
  } else { // return t1 > t2
    return t1 - t2 > std::numeric_limits<T>::epsilon();
  }
}

/*! \return t1 <= t2 if dt > 0 and t1 => t2 if dt < 0 with epsilon accuracy
 *  \note This function is part of Boost.Odeint
 */
template <typename T> bool less_eq_with_sign(T t1, T t2, T dt) {
  if (dt > 0) {
    return t1 - t2 <= std::numeric_limits<T>::epsilon();
  } else {
    return t2 - t1 <= std::numeric_limits<T>::epsilon();
  }
}
} // namespace detail

/*! \brief Integrates an ODE given a stepper.
 *  \author Roberto Di Remigio
 *  \date 2018
 *  \tparam Stepper  integrate stepping function
 *  \tparam System   the ODE system to integrate
 *  \tparam Observer function reporting integration progress
 *
 *  The interface mimics the one of the Boost.Odeint package.
 */
template <typename Stepper, typename System, typename Observer>
size_t integrate_const(const Stepper & stepper,
                       const System & system,
                       typename System::StateType & f,
                       const double tStart,
                       const double tEnd,
                       const double dt,
                       const Observer & observer) {
  size_t nSteps = 0;
  double t = tStart;
  while (detail::less_eq_with_sign(t + dt, tEnd, dt)) {
    observer(f, t);
    stepper.doStep(system, f, t, dt);
    ++nSteps;
    t = tStart + nSteps * dt;
  }
  observer(f, t);
  return nSteps;
}

/*! \brief 4th-order Runge-Kutta stepper.
 *  \author Roberto Di Remigio
 *  \date 2018
 *  \tparam StateType the data structure used to hold the state of the ODE
 *
 *  The interface mimics the one of the Boost.Odeint package.
 */
template <typename StateType> struct RungeKutta4 __final {
  /*! \brief Implements the "classic" 4th-order Runge-Kutta stepper.
   *  \author Roberto Di Remigio
   *  \date 2018
   *  \tparam System the ODE system to integrate
   *
   *  This function implements the "classic" Runge-Kutta stepper, with the
   *  following Butcher Tableau:
   *
   *  \f[
   *    \begin{array}{c|cccc}
   *     0   & 0   & 0   & 0   & 0\\
   *     1/2 & 1/2 & 0   & 0   & 0\\
   *     1/2 & 0   & 1/2 & 0   & 0\\
   *     1   & 0   & 0   & 1   & 0\\
   *     \hline
   *         & 1/6 & 1/3 & 1/3 & 1/6\\
   *    \end{array}
   *  \f]
   */
  template <typename System>
  void doStep(const System & system,
              StateType & f,
              const double t,
              const double dt) const {
    StateType k1, k2, k3, k4;
    StateType x1, x2, x3;

    // Step 1
    system(f, k1, t);
    // Step2
    for (size_t i = 0; i < system.ODEorder(); ++i) {
      x1[i] = f[i] + 0.5 * dt * k1[i];
    }
    system(x1, k2, t + 0.5 * dt);
    // Step 3
    for (size_t i = 0; i < system.ODEorder(); ++i) {
      x2[i] = f[i] + 0.5 * dt * k2[i];
    }
    system(x2, k3, t + 0.5 * dt);
    // Step 4
    for (size_t i = 0; i < system.ODEorder(); ++i) {
      x3[i] = f[i] + dt * k3[i];
    }
    system(x3, k4, t + dt);
    // Update
    for (size_t i = 0; i < system.ODEorder(); ++i) {
      f[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
  }

  /*! \brief Implements the 3/8 4th-order Runge-Kutta stepper.
   *  \author Roberto Di Remigio
   *  \date 2018
   *  \tparam System the ODE system to integrate
   *
   *  This function implements the 3/8 Runge-Kutta stepper, with the
   *  following Butcher Tableau:
   *
   *  \f[
   *    \begin{array}{c|cccc}
   *    0   & 0   & 0   & 0   & 0\\
   *    1/3 & 1/3 & 0   & 0   & 0\\
   *    2/3 & -1/3   & 1 & 0   & 0\\
   *    1   & 1   & -1   & 1   & 0\\
   *    \hline
   *        & 1/8 & 3/8 & 3/8 & 1/8\\
   *    \end{array}
   *  \f]
   */
  template <typename System>
  void doStep38(const System & system,
                StateType & f,
                const double t,
                const double dt) const {
    StateType k1{}, k2{}, k3{}, k4{};
    StateType x1{}, x2{}, x3{};

    // Step 1
    system(f, k1, t);
    // Step 2
    for (size_t i = 0; i < system.ODEorder(); ++i) {
      x1[i] = f[i] + (dt / 3.0) * k1[i];
    }
    system(x1, k2, t + (dt / 3.0));
    // Step 3
    for (size_t i = 0; i < system.ODEorder(); ++i) {
      x2[i] = f[i] + dt * (-1.0 / 3.0 * k1[i] + k2[i]);
    }
    system(x2, k3, t + dt * (2.0 / 3.0));
    // Step 4
    for (size_t i = 0; i < system.ODEorder(); ++i) {
      x3[i] = f[i] + dt * (k1[i] - k2[i] + k3[i]);
    }
    system(x3, k4, t + dt);
    // Update
    for (size_t i = 0; i < system.ODEorder(); ++i) {
      f[i] += dt * ((1.0 / 8.0) * k1[i] + (3.0 / 8.0) * k2[i] + (3.0 / 8.0) * k3[i] +
                    (1.0 / 8.0) * k4[i]);
    }
  }
};
} // namespace utils
} // namespace pcm
