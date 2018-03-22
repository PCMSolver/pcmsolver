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

#include <cmath>
#include <fstream>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#ifndef HAS_CXX11
#include <boost/foreach.hpp>
#endif
// Boost.Odeint includes
#include <boost/numeric/odeint.hpp>

#include "utils/MathUtils.hpp"
#include "utils/RungeKutta4.hpp"

/*! \file InterfacesImpl.hpp */

namespace pcm {
namespace green {
namespace detail {

/*! \typedef ProfileEvaluator
 *  \brief sort of a function pointer to the dielectric profile evaluation function
 */
typedef pcm::function<pcm::tuple<double, double>(const double)> ProfileEvaluator;

/*! \struct IntegratorParameters
 *  \brief holds parameters for the integrator
 */
struct IntegratorParameters {
  /*! Lower bound of the integration interval */
  double r_0_;
  /*! Upper bound of the integration interval */
  double r_infinity_;
  /*! Time step between observer calls */
  double observer_step_;
  IntegratorParameters(double r0, double rinf, double step)
      : r_0_(r0), r_infinity_(rinf), observer_step_(step) {}
};

/*! \class LnTransformedRadial
 *  \brief system of ln-transformed first-order radial differential equations
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  Provides a handle to the system of differential equations for the integrator.
 *  The dielectric profile comes in as a boost::function object.
 */
class LnTransformedRadial __final : public pcm::utils::detail::ODESystem<2> {
public:
  /*! Type of the state vector of the ODE */
  typedef pcm::utils::detail::ODESystem<2>::StateType StateType;
  /*! Constructor from profile evaluator and angular momentum */
  LnTransformedRadial(const ProfileEvaluator & e, int lval) : eval_(e), l_(lval) {}

private:
  /*! Dielectric profile function and derivative evaluation */
  ProfileEvaluator eval_;
  /*! Angular momentum */
  int l_;
  /*! Provides a functor for the evaluation of the system
   *  of first-order ODEs needed by Boost.Odeint
   *  The second-order ODE and the system of first-order ODEs
   *  are reported in the manuscript.
   *  \param[in] rho state vector holding the function and its first derivative
   *  \param[out] drhodr state vector holding the first and second derivative
   *  \param[in] y logarithmic position on the integration grid
   */
  virtual void RHS(const StateType & rho, StateType & drhodr, const double y) const {
    // Evaluate the dielectric profile
    double eps = 0.0, epsPrime = 0.0;
    pcm::tie(eps, epsPrime) = eval_(std::exp(y));
    if (utils::numericalZero(eps))
      PCMSOLVER_ERROR("Division by zero!");
    double gamma_epsilon = std::exp(y) * epsPrime / eps;
    // System of equations is defined here
    drhodr[0] = rho[1];
    drhodr[1] = -rho[1] * (rho[1] + 1.0 + gamma_epsilon) + l_ * (l_ + 1);
  }
};
} // namespace detail

using detail::IntegratorParameters;
using detail::ProfileEvaluator;

/*! \class RadialFunction
 *  \brief represents solutions to the radial 2nd order ODE
 *  \author Roberto Di Remigio
 *  \date 2018
 *  \tparam ODE system of 1st order ODEs replacing the 2nd order ODE
 *
 * The sign of the integration interval, i.e. the sign of the difference
 * between the minimum and the maximum of the interval, is used to discriminate
 * between the two independent solutions (zeta-type and omega-type) of the
 * differential equation:
 *   - When the sign is positive (y_min_ < y_max_), we are integrating **forwards**
 * for the zeta-type solution.
 *   - When the sign is negative (y_max_ < y_min_), we are integrating **backwards**
 * for the omega-type solution.
 *
 * This is based on an idea Luca had in Wuhan/Chongqing for the planar
 * diffuse interface.
 * With respect to the previous, much more heavily template-d implementation,
 * it introduces a bit more if-else if-else branching in the function_impl and
 * derivative_impl functions.
 * This is largely offset by the better readability and reduced code duplication.
 */
template <typename ODE> class RadialFunction __final {
public:
  RadialFunction() : L_(0), y_min_(0.0), y_max_(0.0), y_sign_(1) {}
  RadialFunction(int l,
                 double ymin,
                 double ymax,
                 const ProfileEvaluator & eval,
                 const IntegratorParameters & parms)
      : L_(l),
        y_min_(ymin),
        y_max_(ymax),
        y_sign_(pcm::utils::sign(y_max_ - y_min_)) {
    compute(eval, parms);
  }
  pcm::tuple<double, double> operator()(double point) const {
    return pcm::make_tuple(function_impl(point), derivative_impl(point));
  }
  friend std::ostream & operator<<(std::ostream & os, RadialFunction & obj) {
    for (size_t i = 0; i < obj.function_[0].size(); ++i) {
      os << obj.function_[0][i] << "    " << obj.function_[1][i] << "    "
         << obj.function_[2][i] << std::endl;
    }
    return os;
  }

private:
  typedef typename ODE::StateType StateType;
  /*! Angular momentum of the solution */
  int L_;
  /*! Lower bound of the integration interval */
  double y_min_;
  /*! Upper bound of the integration interval */
  double y_max_;
  /*! The sign of the integration interval */
  double y_sign_;
  /*! The actual data: grid, function value and first derivative values */
  pcm::array<std::vector<double>, 3> function_;
  /*! Reports progress of differential equation integrator */
  void push_back(const StateType & x, double y) {
    function_[0].push_back(y);
    function_[1].push_back(x[0]);
    function_[2].push_back(x[1]);
  }
  /*! \brief Calculates radial solution
   *  \param[in] eval   dielectric profile evaluator function object
   *  \param[in] parms parameters for the integrator
   *
   *  This function discriminates between the first, i.e. the one with r^l
   *  behavior, and the second radial solution, i.e. the one with r^(-l-1)
   *  behavior, based on the sign of the integration interval y_sign_.
   */
  void compute(const ProfileEvaluator & eval, const IntegratorParameters & parms) {
    namespace odeint = boost::numeric::odeint;
    odeint::runge_kutta4<StateType> stepper;

    ODE system(eval, L_);
    // Holds the initial conditions
    StateType init;
    // Set initial conditions
    if (y_sign_ > 0.0) { // zeta-type solution
      init[0] = y_sign_ * L_ * y_min_;
      init[1] = y_sign_ * L_;
    } else { // omega-type solution
      init[0] = y_sign_ * (L_ + 1) * y_min_;
      init[1] = y_sign_ * (L_ + 1);
    }
    odeint::integrate_const(
        stepper,
        system,
        init,
        y_min_,
        y_max_,
        y_sign_ * parms.observer_step_,
        pcm::bind(&RadialFunction<ODE>::push_back, this, pcm::_1, pcm::_2));
    // clang-format off
    // Reverse order of function_ if omega-type solution was computed
    // this ensures that they are in ascending order, as later expected by
    // function_impl and derivative_impl
    if (y_sign_ < 0.0) {
#ifdef HAS_CXX11
      for (auto & comp : function_) {
#else  /* HAS_CXX11 */
      BOOST_FOREACH (std::vector<double> & comp, function_) {
#endif /* HAS_CXX11 */
        std::reverse(comp.begin(), comp.end());
      }
      }
    }
  // clang-format on
  /*! \brief Returns value of function at given point
   *  \param[in] point evaluation point
   *
   *  We first check if point is below y_min_, if yes we use
   *  the asymptotic form.
   */
  double function_impl(double point) const {
    double val = 0.0;
    if (point <= y_min_ && y_sign_ > 0.0) { // Asymptotic zeta-type solution
      val = y_sign_ * L_ * point;
    } else if (point >= y_min_ && y_sign_ < 0.0) { // Asymptotic omega-type solution
      val = y_sign_ * (L_ + 1) * point;
    } else {
      val = utils::splineInterpolation(point, function_[0], function_[1]);
    }
    return val;
  }
  /*! \brief Returns value of 1st derivative of function at given point
   *  \param[in] point evaluation point
   *
   *  Below y_min_, the asymptotic form is used. Otherwise we interpolate.
   */
  double derivative_impl(double point) const {
    double val = 0.0;
    if (point <= y_min_ && y_sign_ > 0.0) { // Asymptotic zeta-type solution
      val = y_sign_ * L_;
    } else if (point >= y_min_ && y_sign_ < 0.0) { // Asymptotic omega-type solution
      val = y_sign_ * (L_ + 1);
    } else {
      val = utils::splineInterpolation(point, function_[0], function_[2]);
    }
    return val;
  }
};

/*! \brief Write contents of a RadialFunction to file
 *  \param[in] f solution to the radial 2nd order ODE to print
 *  \param[in] fname name of the file
 *  \author Roberto Di Remigio
 *  \date 2018
 *  \tparam ODE system of 1st order ODEs replacing the 2nd order ODE
 */
template <typename ODE>
void writeToFile(RadialFunction<ODE> & f, const std::string & fname) {
  std::ofstream fout;
  fout.open(fname.c_str());
  fout << f << std::endl;
  fout.close();
}
} // namespace green
} // namespace pcm
