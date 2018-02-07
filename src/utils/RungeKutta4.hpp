/*
 * The Runge-Kutta source code is Copyright(c) 2013 John Burkardt.
 *
 * This Runge-Kutta ODE solver is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This Runge-Kutts ODE Solver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the code.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See also:
 * https://people.sc.fsu.edu/~jburkardt/cpp_src/rk4/rk4.html
 */

/*! \file RungeKutta4.hpp */

namespace pcm {
namespace utils {
/*! Right-hand side for a scalar ordinary-differential equation */
typedef std::function<double(double, double)> ScalarODERHS;

/*! \brief  RK4 takes one Runge-Kutta step for a scalar ODE.
 *  \author John Burkardt, Roberto Di Remigio
 *  \param[in] t0 the current time.
 *  \param[in] u0 the solution estimate at the current time.
 *  \param[in] dt the time step.
 *  \param[in] F a function which evaluates the derivative, or right hand side of the
 * problem.
 *  \return the fourth-order Runge-Kutta solution estimate at time t0+dt.
 *
 *  It is assumed that an initial value problem, of the form
 *
 *    du/dt = f( t, u )
 *    u(t0) = u0
 *
 *  is being solved.
 *
 *  If the user can supply current values of t, u, a stepsize dt, and a
 *  function to evaluate the derivative, this function can compute the
 *  fourth-order Runge Kutta estimate to the solution at time t+dt.
 */
double rk4(double t0, double u0, double dt, double f(double t, double u));

/*! Right-hand side for a vector ordinary-differential equation */
typedef std::function<double *(double, int, double)> VectorODERHS;

/*!
 *  \brief Takes one Runge-Kutta step for a vector ODE.
 *  \author John Burkardt, Roberto Di Remigio
 *  \param[in] t0 the current time.
 *  \param[in] n the spatial dimension.
 *  \param[in] u0[M] the solution estimate at the current time.
 *  \param[in] dt the time step.
 *  \param[in] F a function which evaluates the derivative, or right hand side of the
 * problem.
 *  \return the fourth-order Runge-Kutta solution estimate at time t0+dt.
 *
 *  It is assumed that an initial value problem, of the form
 *
 *    du/dt = f ( t, u )
 *    u(t0) = u0
 *
 *  is being solved.
 *
 *  If the user can supply current values of t, u, a stepsize dt, and a
 *  function to evaluate the derivative, this function can compute the
 *  fourth-order Runge Kutta estimate to the solution at time t+dt.
 */
double * rk4vec(double t0,
                int n,
                double u0[],
                double dt,
                double * f(double t, int n, double u[]));
} // namespace utils
} // namespace pcm
