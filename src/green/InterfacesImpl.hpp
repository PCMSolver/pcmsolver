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

#include <cmath>
#include <fstream>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>
#include <boost/foreach.hpp>

// Boost.Odeint includes
#include <boost/numeric/odeint.hpp>

#include "utils/MathUtils.hpp"

namespace pcm {
namespace green {
namespace detail {

/*! \typedef StateType
 *  \brief state vector for the differential equation integrator
 */
typedef std::vector<double> StateType;

/*! \typedef ProfileEvaluator
 *  \brief sort of a function pointer to the dielectric profile evaluation function
 */
typedef pcm::function<pcm::tuple<double, double>(const double)> ProfileEvaluator;

/*! \struct IntegratorParameters
 *  \brief holds parameters for the integrator
 */
struct IntegratorParameters {
  /*! Absolute tolerance level */
  double eps_abs_;
  /*! Relative tolerance level */
  double eps_rel_;
  /*! Weight of the state      */
  double factor_x_;
  /*! Weight of the state derivative */
  double factor_dxdt_;
  /*! Lower bound of the integration interval */
  double r_0_;
  /*! Upper bound of the integration interval */
  double r_infinity_;
  /*! Time step between observer calls */
  double observer_step_;
  IntegratorParameters(double e_abs,
                       double e_rel,
                       double f_x,
                       double f_dxdt,
                       double r0,
                       double rinf,
                       double step)
      : eps_abs_(e_abs),
        eps_rel_(e_rel),
        factor_x_(f_x),
        factor_dxdt_(f_dxdt),
        r_0_(r0),
        r_infinity_(rinf),
        observer_step_(step) {}
};

/*! \class LnTransformedRadial
 *  \brief system of ln-transformed first-order radial differential equations
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  Provides a handle to the system of differential equations for the integrator.
 *  The dielectric profile comes in as a boost::function object.
 */
class LnTransformedRadial {
private:
  /*! Dielectric profile function and derivative evaluation */
  ProfileEvaluator eval_;
  /*! Angular momentum */
  int l_;

public:
  /*! Constructor from profile evaluator and angular momentum */
  LnTransformedRadial(const ProfileEvaluator & e, int lval) : eval_(e), l_(lval) {}
  /*! Provides a functor for the evaluation of the system
   *  of first-order ODEs needed by Boost.Odeint
   *  The second-order ODE and the system of first-order ODEs
   *  are reported in the manuscript.
   *  \param[in] rho state vector holding the function and its first derivative
   *  \param[out] drhodr state vector holding the first and second derivative
   *  \param[in] r position on the integration grid
   */
  void operator()(const StateType & rho, StateType & drhodr, const double r) {
    // Evaluate the dielectric profile
    double eps = 0.0, epsPrime = 0.0;
    pcm::tie(eps, epsPrime) = eval_(r);
    if (utils::numericalZero(eps))
      throw std::domain_error("Division by zero!");
    double gamma_epsilon = epsPrime / eps;
    // System of equations is defined here
    drhodr[0] = rho[1];
    drhodr[1] = -rho[1] * (rho[1] + 2.0 / r + gamma_epsilon) +
                l_ * (l_ + 1) / std::pow(r, 2);
  }
};
} // namespace detail

using detail::ProfileEvaluator;
using detail::IntegratorParameters;

/*! \file InterfacesImpl.hpp
 *  \class RadialFunction
 *  \brief represents solutions to the radial 2nd order ODE
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam StateVariable type of the state variable used in the ODE solver
 *  \tparam ODESystem system of 1st order ODEs replacing the 2nd order ODE
 *  \tparam IndependentSolution encodes which type of radial solution
 */
template <typename StateVariable,
          typename ODESystem,
          template <typename, typename> class IndependentSolution>
class RadialFunction __final {
public:
  RadialFunction() : solution_(IndependentSolution<StateVariable, ODESystem>()) {}
  RadialFunction(int l,
                 double r0,
                 double rinf,
                 const ProfileEvaluator & eval,
                 const IntegratorParameters & parms)
      : solution_(IndependentSolution<StateVariable, ODESystem>(l,
                                                                r0,
                                                                rinf,
                                                                eval,
                                                                parms)) {}
  ~RadialFunction() {}
  /*! \brief Returns value of function and its first derivative at given point
   *  \param[in] point evaluation point
   */
  pcm::tuple<double, double> operator()(double point) const {
    return solution_(point);
  }
  friend std::ostream & operator<<(std::ostream & os, RadialFunction & obj) {
    os << obj.solution_;
    return os;
  }

private:
  /// Independent solution to the radial equation
  IndependentSolution<StateVariable, ODESystem> solution_;
};

/*! \file InterfacesImpl.hpp
 *  \class Zeta
 *  \brief 1st solution to the radial second order ODE, with r^l behaviour in the
 * origin
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam StateVariable type of the state variable used in the ODE solver
 *  \tparam ODESystem system of 1st order ODEs replacing the 2nd order ODE
 */
template <typename StateVariable, typename ODESystem> class Zeta __final {
public:
  Zeta() : L_(0), r_0_(0.0), r_infinity_(0.0) {}
  Zeta(int l,
       double r0,
       double rinf,
       const ProfileEvaluator & eval,
       const IntegratorParameters & parms)
      : L_(l), r_0_(r0), r_infinity_(rinf) {
    compute(eval, parms);
  }
  ~Zeta() {}
  pcm::tuple<double, double> operator()(double point) const {
    return pcm::make_tuple(function_impl(point), derivative_impl(point));
  }
  friend std::ostream & operator<<(std::ostream & os, Zeta & obj) {
    for (size_t i = 0; i < obj.function_[0].size(); ++i) {
      os << obj.function_[0][i] << "    " << obj.function_[1][i] << "    "
         << obj.function_[2][i] << std::endl;
    }
    return os;
  }

private:
  typedef pcm::array<StateVariable, 3> RadialSolution;
  /// Angular momentum of the solution
  int L_;
  /// Lower bound of the integration interval
  double r_0_;
  /// Upper bound of the integration interval
  double r_infinity_;
  /// The actual data: grid, function value and first derivative values
  RadialSolution function_;
  /*! Reports progress of differential equation integrator */
  void push_back(const StateVariable & x, double r) {
    function_[0].push_back(r);
    function_[1].push_back(x[0]);
    function_[2].push_back(x[1]);
  }
  /*! \brief Calculates 1st radial solution, i.e. the one with r^l behavior
   *  \param[in] eval   dielectric profile evaluator function object
   *  \param[in] parms parameters for the integrator
   */
  void compute(const ProfileEvaluator & eval, const IntegratorParameters & parms) {
    namespace odeint = boost::numeric::odeint;
    odeint::bulirsch_stoer_dense_out<StateVariable> stepper(
        parms.eps_abs_, parms.eps_rel_, parms.factor_x_, parms.factor_dxdt_);
    ODESystem system(eval, L_);
    // Holds the initial conditions
    StateVariable init_zeta(2);
    // Set initial conditions
    init_zeta[0] = L_ * std::log(r_0_);
    init_zeta[1] = L_ / r_0_;
    odeint::integrate_adaptive(
        stepper,
        system,
        init_zeta,
        r_0_,
        r_infinity_,
        parms.observer_step_,
        pcm::bind(
            &Zeta<StateVariable, ODESystem>::push_back, this, pcm::_1, pcm::_2));
  }
  /*! \brief Returns value of function at given point
   *  \param[in] point evaluation point
   *
   *  We first check if point is below r_0_, if yes we use
   *  the asymptotic form L*log(r) in point.
   */
  double function_impl(double point) const {
    double zeta = 0.0;
    if (point <= r_0_) {
      zeta = L_ * std::log(point);
    } else {
      zeta = utils::splineInterpolation(point, function_[0], function_[1]);
    }
    return zeta;
  }
  /*! \brief Returns value of 1st derivative of function at given point
   *  \param[in] point evaluation point
   *
   *  We first check if point is below r_0_, if yes we use
   *  the asymptotic form L / r in point.
   */
  double derivative_impl(double point) const {
    double zeta = 0.0;
    if (point <= r_0_) {
      zeta = L_ / point;
    } else {
      zeta = utils::splineInterpolation(point, function_[0], function_[2]);
    }
    return zeta;
  }
};

/*! \file InterfacesImpl.hpp
 *  \class Omega
 *  \brief 2nd solution to the radial second order ODE, with r^(-l-1) behaviour at
 * infinity
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam StateVariable type of the state variable used in the ODE solver
 *  \tparam ODESystem system of 1st order ODEs replacing the 2nd order ODE
 */
template <typename StateVariable, typename ODESystem> class Omega __final {
public:
  Omega() : L_(0), r_0_(0.0), r_infinity_(0.0) {}
  Omega(int l,
        double r0,
        double rinf,
        const ProfileEvaluator & eval,
        const IntegratorParameters & parms)
      : L_(l), r_0_(r0), r_infinity_(rinf) {
    compute(eval, parms);
  }
  ~Omega() {}
  pcm::tuple<double, double> operator()(double point) const {
    return pcm::make_tuple(function_impl(point), derivative_impl(point));
  }
  friend std::ostream & operator<<(std::ostream & os, Omega & obj) {
    for (size_t i = 0; i < obj.function_[0].size(); ++i) {
      os << obj.function_[0][i] << "    " << obj.function_[1][i] << "    "
         << obj.function_[2][i] << std::endl;
    }
    return os;
  }

private:
  typedef pcm::array<StateVariable, 3> RadialSolution;
  /// Angular momentum of the solution
  int L_;
  /// Lower bound of the integration interval
  double r_0_;
  /// Upper bound of the integration interval
  double r_infinity_;
  /// The actual data: grid, function value and first derivative values
  RadialSolution function_;
  /*! Reports progress of differential equation integrator */
  void push_back(const StateVariable & x, double r) {
    function_[0].push_back(r);
    function_[1].push_back(x[0]);
    function_[2].push_back(x[1]);
  }
  /*! \brief calculates 2nd radial solution, i.e. the one with r^(-l-1) behavior
   *  \param[in] eval   dielectric profile evaluator function object
   *  \param[in] parms parameters for the integrator
   */
  void compute(const ProfileEvaluator & eval, const IntegratorParameters & parms) {
    namespace odeint = boost::numeric::odeint;
    odeint::bulirsch_stoer_dense_out<StateVariable> stepper(
        parms.eps_abs_, parms.eps_rel_, parms.factor_x_, parms.factor_dxdt_);
    ODESystem system(eval, L_);
    // Holds the initial conditions
    StateVariable init_omega(2);
    // Set initial conditions
    init_omega[0] = -(L_ + 1) * std::log(r_infinity_);
    init_omega[1] = -(L_ + 1) / r_infinity_;
    // Notice that we integrate BACKWARDS, so we pass -step to integrate_adaptive
    boost::numeric::odeint::integrate_adaptive(
        stepper,
        system,
        init_omega,
        r_infinity_,
        r_0_,
        -parms.observer_step_,
        pcm::bind(
            &Omega<StateVariable, ODESystem>::push_back, this, pcm::_1, pcm::_2));
    // Reverse order of StateVariable-s in RadialSolution
    // this ensures that they are in ascending order, as later expected by
    // function_impl and derivative_impl
    BOOST_FOREACH (StateVariable & comp, function_) {
      std::reverse(comp.begin(), comp.end());
    }
  }
  /*! \brief Returns value of function at given point
   *  \param[in] point evaluation point
   *
   * We first check if point is above r_infinity_, if yes we use
   * the asymptotic form -(L+1)*log(r) in point.
   */
  double function_impl(double point) const {
    double omega = 0.0;
    if (point >= r_infinity_) {
      omega = -(L_ + 1) * std::log(point);
    } else {
      omega = utils::splineInterpolation(point, function_[0], function_[1]);
    }
    return omega;
  }
  /*! \brief Returns value of 1st derivative of function at given point
   *  \param[in] point evaluation point
   *
   * We first check if point is above r_infinity_, if yes we use
   * the asymptotic form -(L+1)/r in point.
   */
  double derivative_impl(double point) const {
    double omega = 0.0;
    if (point >= r_infinity_) {
      omega = -(L_ + 1) / point;
    } else {
      omega = utils::splineInterpolation(point, function_[0], function_[2]);
    }
    return omega;
  }
};

/*! \brief Write contents of a RadialFunction to file
 *  \param[in] f RadialSolution whose contents have to be printed
 *  \param[in] fname name of the file
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam StateVariable type of the state variable used in the ODE solver
 *  \tparam ODESystem system of 1st order ODEs replacing the 2nd order ODE
 *  \tparam IndependentSolution encodes which type of radial solution
 */
template <typename StateVariable,
          typename ODESystem,
          template <typename, typename> class IndependentSolution>
void writeToFile(RadialFunction<StateVariable, ODESystem, IndependentSolution> & f,
                 const std::string & fname) {
  std::ofstream fout;
  fout.open(fname.c_str());
  fout << f << std::endl;
  fout.close();
}
} // namespace green
} // namespace pcm
