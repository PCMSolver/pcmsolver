/*
 * The Runge-Kutta source code is Copyright(c) 2010 John Burkardt.
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

#include "rk4.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace pcm {
namespace utils {
using namespace std;

double rk4(double t0, double u0, double dt, double f(double t, double u)) {
  double f0;
  double f1;
  double f2;
  double f3;
  double t1;
  double t2;
  double t3;
  double u;
  double u1;
  double u2;
  double u3;
  //
  //  Get four sample values of the derivative.
  //
  f0 = f(t0, u0);

  t1 = t0 + dt / 2.0;
  u1 = u0 + dt * f0 / 2.0;
  f1 = f(t1, u1);

  t2 = t0 + dt / 2.0;
  u2 = u0 + dt * f1 / 2.0;
  f2 = f(t2, u2);

  t3 = t0 + dt;
  u3 = u0 + dt * f2;
  f3 = f(t3, u3);
  //
  //  Combine to estimate the solution at time T0 + DT.
  //
  u = u0 + dt * (f0 + 2.0 * f1 + 2.0 * f2 + f3) / 6.0;

  return u;
}

double * rk4vec(double t0,
                int m,
                double u0[],
                double dt,
                double * f(double t, int m, double u[]))

{
  double * f0;
  double * f1;
  double * f2;
  double * f3;
  int i;
  double t1;
  double t2;
  double t3;
  double * u;
  double * u1;
  double * u2;
  double * u3;
  //
  //  Get four sample values of the derivative.
  //
  f0 = f(t0, m, u0);

  t1 = t0 + dt / 2.0;
  u1 = new double[m];
  for (i = 0; i < m; i++) {
    u1[i] = u0[i] + dt * f0[i] / 2.0;
  }
  f1 = f(t1, m, u1);

  t2 = t0 + dt / 2.0;
  u2 = new double[m];
  for (i = 0; i < m; i++) {
    u2[i] = u0[i] + dt * f1[i] / 2.0;
  }
  f2 = f(t2, m, u2);

  t3 = t0 + dt;
  u3 = new double[m];
  for (i = 0; i < m; i++) {
    u3[i] = u0[i] + dt * f2[i];
  }
  f3 = f(t3, m, u3);
  //
  //  Combine them to estimate the solution.
  //
  u = new double[m];
  for (i = 0; i < m; i++) {
    u[i] = u0[i] + dt * (f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i]) / 6.0;
  }
  //
  //  Free memory.
  //
  delete[] f0;
  delete[] f1;
  delete[] f2;
  delete[] f3;
  delete[] u1;
  delete[] u2;
  delete[] u3;

  return u;
}
} // namespace utils
} // namespace pcm
