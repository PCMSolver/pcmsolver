/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef TANHSPHERICALDIFFUSE_HPP
#define TANHSPHERICALDIFFUSE_HPP

#include <cmath>
#include <functional>
#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>
#include "taylor.hpp"

// Boost.Odeint includes
#include <boost/numeric/odeint/config.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>
// Boost.Math includes
#include <boost/math/special_functions/legendre.hpp>

#include "DerivativeTypes.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"
#include "SphericalDiffuse.hpp"
#include "TanhDiffuse.hpp"
#include "Timer.hpp"
#include "LoggerInterface.hpp"
#include "MathUtils.hpp"

/*! \file TanhSphericalDiffuse.hpp
 *  \class LnTransformedRadial
 *  \brief system of ln-transformed first-order radial differential equations
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  Provides a handle to the system of differential equations for the integrator.
 *  The dielectric profile comes in as a boost::function object.
 */
class LnTransformedRadial
{
    private:
        /*! Dielectric profile function and derivative evaluation */
        ProfileEvaluator eval_;
        /*! Angular momentum */
        size_t l_;
    public:
        /*! Constructor from profile evaluator and angular momentum */
        LnTransformedRadial(const ProfileEvaluator & e, size_t lval) : eval_(e), l_(lval) {}
        /*! Provides a functor for the evaluation of the system
         *  of first-order ODEs needed by Boost.Odeint
         *  The second-order ODE and the system of first-order ODEs
         *  are reported in the manuscript.
         *  \param[in] rho state vector holding the function and its first derivative
         *  \param[out] drhodr state vector holding the first and second derivative
         *  \param[in] r position on the integration grid
         */
        void operator()(const StateType & rho, StateType & drhodr, const double r)
        {
            // Evaluate the dielectric profile
            double eps = 0.0, epsPrime = 0.0;
            eval_(eps, epsPrime, r);
            double gamma_epsilon = epsPrime / eps;
            // System of equations is defined here
            drhodr[0] = rho[1];
            drhodr[1] = -rho[1] * (rho[1] + 2.0/r + gamma_epsilon) + l_ * (l_ + 1) / std::pow(r, 2) * rho[0];
        }
};

/*! \file TanhSphericalDiffuse.hpp
 *  \brief reports progress of differential equation integrator
 *  \author Roberto Di Remigio
 *  \date 2015
 */
inline void observer(RadialFunction & f, const StateType & x, double r)
{
    /* Save grid points */
    f[0].push_back(r);
    /* Save function */
    f[1].push_back(x[0]);
    /* Save first derivative of function */
    f[2].push_back(x[1]);
}

/*! \brief Calculates 1st radial solution, i.e. the one with r^l behavior
 *  \param[in]  L      angular momentum of the required solution
 *  \param[out] f      solution to the radial equation
 *  \param[in]  eval   dielectric profile evaluator function object
 *  \param[in]  params parameters for the integrator
 */
inline void computeZeta(const size_t L, RadialFunction & f, const ProfileEvaluator & eval, const IntegratorParameters & params)
{
    using namespace std::placeholders;
    namespace odeint = boost::numeric::odeint;
    odeint::bulirsch_stoer_dense_out<StateType> stepper_(params.eps_abs_, params.eps_rel_, params.factor_x_, params.factor_dxdt_);
    // The system of first-order ODEs
    LnTransformedRadial system_(eval, L);
    // Holds the initial conditions
    StateType init_zeta_(2);
    // Partial application of observer
    auto observer_ = std::bind(observer, f, _1, _2);
    // Set initial conditions
    init_zeta_[0] = L * std::log(params.r_0_);
    init_zeta_[1] = L / params.r_0_;
    odeint::integrate_adaptive(stepper_, system_, init_zeta_, params.r_0_, params.r_infinity_, params.observer_step_, observer_);
}

/*! \brief Calculates 2nd radial solution, i.e. the one with r^(-l-1) behavior
 *  \param[in]  L      angular momentum of the required solution
 *  \param[out] f      solution to the radial equation
 *  \param[in]  eval   dielectric profile evaluator function object
 *  \param[in]  params parameters for the integrator
 */
inline void computeOmega(const size_t L, RadialFunction & f, const ProfileEvaluator & eval, const IntegratorParameters & params)
{
    using namespace std::placeholders;
    namespace odeint = boost::numeric::odeint;
    odeint::bulirsch_stoer_dense_out<StateType> stepper_(params.eps_abs_, params.eps_rel_, params.factor_x_, params.factor_dxdt_);
    // The system of first-order ODEs
    LnTransformedRadial system_(eval, L);
    // Holds the initial conditions
    StateType init_omega_(2);
    // Partial application of observer
    auto observer_ = std::bind(observer, f, _1, _2);
    // Set initial conditions
    init_omega_[0] = -(L + 1) * std::log(params.r_infinity_);
    init_omega_[1] = -(L + 1) / params.r_infinity_;
    // Notice that we integrate BACKWARDS, so we pass -params.observer_step_ to integrate_adaptive
    odeint::integrate_adaptive(stepper_, system_, init_omega_, params.r_infinity_, params.r_0_, -params.observer_step_, observer_);
}

template <>
inline void TanhSphericalDiffuse::initSphericalDiffuse()
{
    using namespace std::placeholders;

    LOG("TanhSphericalDiffuse::initSphericalDiffuse");
    // Parameters for the numerical solution of the radial differential equation
    // Initialize the system of differential equations
    double eps_abs_     = 1.0e-10; /*! Absolute tolerance level */
    double eps_rel_     = 1.0e-06; /*! Relative tolerance level */
    double factor_x_    = 0.0;     /*! Weight of the state      */
    double factor_dxdt_ = 0.0;     /*! Weight of the state derivative */
    double r_0_         = 0.5;     /*! Lower bound of the integration interval */
    double r_infinity_  = profile_.center() + 100.0; /*! Upper bound of the integration interval */
    double observer_step_ = 5.0e-04; /*! Time step between observer calls */
    IntegratorParameters params_(eps_abs_, eps_rel_, factor_x_, factor_dxdt_, r_0_, r_infinity_, observer_step_);
    ProfileEvaluator eval_ = std::bind(&TanhDiffuse::operator(), this->profile_, _1, _2, _3);

    LOG("Computing coefficient for the separation of the Coulomb singularity");
    LOG("Computing first radial solution L = ", maxLC_);
    timerON("TanhSphericalDiffuse::computeZeta for coefficient");
    computeZeta(maxLC_, zetaC_, eval_, params_);
    timerOFF("TanhSphericalDiffuse::computeZeta for coefficient");
    LOG("DONE: Computing first radial solution L = ", maxLC_);

    LOG("Computing first radial solution L = ", maxLC_);
    timerON("TanhSphericalDiffuse::computeOmega");
    computeOmega(maxLC_, omegaC_, eval_, params_);
    timerOFF("TanhSphericalDiffuse::computeOmega");
    LOG("Computing second radial solution L = ", maxLC_);
    LOG("DONE: Computing coefficient for the separation of the Coulomb singularity");

    LOG("Computing radial solutions for Green's function");
    timerON("   Looping over angular momentum");
    for (size_t L = 0; L < maxLGreen_; ++L) {
        // First radial solution
        LOG("Computing first radial solution L = ", L);
        timerON("TanhSphericalDiffuse::computeZeta");
        // Create an empty RadialFunction
        RadialFunction tmp_zeta_;
        computeZeta(L, tmp_zeta_, eval_, params_);
        zeta_.push_back(tmp_zeta_);
        timerOFF("TanhSphericalDiffuse::computeZeta");
        LOG("DONE: Computing first radial solution L = ", L);

        // Second radial solution
        LOG("Computing first radial solution L = ", L);
        timerON("TanhSphericalDiffuse::computeOmega");
        // Create an empty RadialFunction
        RadialFunction tmp_omega_;
        computeOmega(L, tmp_omega_, eval_, params_);
        omega_.push_back(tmp_omega_);
        timerOFF("TanhSphericalDiffuse::computeOmega");
        LOG("Computing second radial solution L = ", L);
    }
    timerOFF("   Looping over angular momentum");
    LOG("DONE: Computing radial solutions for Green's function");
}

template <>
inline double TanhSphericalDiffuse::coefficient(double r1, double r2) const
{
    /* Value of zetaC_ at point with index 1 */
    double zeta1  = linearInterpolation(r1, zetaC_[0], zetaC_[1]);
    /* Value of zetaC_ at point with index 2 */
    double zeta2  = linearInterpolation(r2, zetaC_[0], zetaC_[1]);
    /* Value of omegaC_ at point with index 1 */
    double omega1 = linearInterpolation(r1, omegaC_[0], omegaC_[1]);
    /* Value of omegaC_ at point with index 2 */
    double omega2 = linearInterpolation(r2, omegaC_[0], omegaC_[1]);

    /* Components for the evaluation of the Wronskian */
    /* Value of derivative of zetaC_ at point with index 2 */
    double d_zeta2  = linearInterpolation(r2, zetaC_[0], zetaC_[2]);
    /* Value of derivative of omegaC_ at point with index 2 */
    double d_omega2 = linearInterpolation(r2, omegaC_[0], omegaC_[2]);

    double tmp = 0.0, coeff = 0.0;
    double eps_r2 = 0.0, epsPrime_r2 = 0.0;
    this->profile_(eps_r2, epsPrime_r2, r2);

    double denominator = (d_zeta2 - d_omega2) * std::pow(r2, 2) * eps_r2;

    if (r1 < r2) {
        tmp = std::exp(zeta1 - zeta2) * (2*maxLC_ +1) / denominator;
        coeff = std::pow(r1/r2, maxLC_) / (tmp * r2);
    } else {
        tmp = std::exp(omega1 - omega2) * (2*maxLC_ +1) / denominator;
        coeff = std::pow(r2/r1, maxLC_) / (tmp * r1);
    }

    return coeff;
}

template <>
inline double TanhSphericalDiffuse::functionSummation(size_t L, double r1, double r2, double cos_gamma, double Cr12) const
{
    double gr12 = 0.0;
    // Evaluate Legendre polynomial of order L
    double pl_x = boost::math::legendre_p(L, cos_gamma);

    /* Value of zeta_[L] at point with index 1 */
    double zeta1  = linearInterpolation(r1, zeta_[L][0], zeta_[L][1]);
    /* Value of zeta_[L} at point with index 2 */
    double zeta2  = linearInterpolation(r2, zeta_[L][0], zeta_[L][1]);
    /* Value of omega_[L] at point with index 1 */
    double omega1 = linearInterpolation(r1, omega_[L][0], omega_[L][1]);
    /* Value of omega_[L} at point with index 2 */
    double omega2 = linearInterpolation(r2, omega_[L][0], omega_[L][1]);

    /* Components for the evaluation of the Wronskian */
    /* Value of derivative of zeta_[L] at point with index 2 */
    double d_zeta2  = linearInterpolation(r2, zeta_[L][0], zeta_[L][2]);
    /* Value of derivative of omega_[L] at point with index 2 */
    double d_omega2 = linearInterpolation(r2, omega_[L][0], omega_[L][2]);

    double eps_r2 = 0.0, epsPrime_r2 = 0.0;
    this->profile_(eps_r2, epsPrime_r2, r2);

    double denominator = (d_zeta2 - d_omega2) * std::pow(r2, 2) * eps_r2;

    if (r1 < r2) {
        gr12 = std::exp(zeta1 - zeta2) * (2*L +1) / denominator;
        gr12 = (gr12 - std::pow(r1/r2, L) / (r2 * Cr12) ) * pl_x ;
    } else {
        gr12 = std::exp(omega1 - omega2) * (2*L +1) / denominator;
        gr12 = (gr12 - std::pow(r2/r1, L) / (r1 * Cr12) ) * pl_x ;
    }

    return gr12;
}

/*! Calcualtes the Green's function given a pair of points */
template <>
inline Numerical GreensFunction<Numerical, TanhDiffuse>::function(const Eigen::Vector3d & source,
                                        const Eigen::Vector3d & probe) const
{
    Numerical sp[3], pp[3];
    sp[0] = source(0); sp[1] = source(1); sp[2] = source(2);
    pp[0] = probe(0);  pp[1] = probe(1);  pp[2] = probe(2);
    return this->operator()(sp, pp);
}


/* Derivative of Legendre polynomials can be obtained in terms of the associated
 * Legendre polynomials (see http://en.wikipedia.org/wiki/Associated_Legendre_polynomials)
 * Condon-Shortley phase is assumed.
 *
 * This is the source code for DOLFIN
 * First derivative of Legendre polynomial of order l at point x:
 * double pl_first_x = -boost::math::legendre_p(l, 1, x) / (std::sqrt(1.0 - x*x));
 * Second derivative of Legendre polynomial of order l at point x:
 * double pl_second_x = boost::math::legendre_p(l, 2, x)/ (1.0 - x*x);
 */

template <>
inline Numerical GreensFunction<Numerical, TanhDiffuse>::derivativeSource(const Eigen::Vector3d & normal_p1,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    // To be implemented
    return 0.0;
}

template <>
inline Numerical GreensFunction<Numerical, TanhDiffuse>::derivativeProbe(const Eigen::Vector3d & normal_p2,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    // To be implemented
    return 0.0;
}

namespace
{
    IGreensFunction * createTanhSphericalDiffuse(const greenData & _data)
    {
        double eL, eR, w, c; // To be read from _data
        eL = 0.0, eR = 0.0, w = 0.0, c = 0.0;
        DiagonalIntegrator * integrator =
		DiagonalIntegratorFactory::TheDiagonalIntegratorFactory().createDiagonalIntegrator(_data.integratorType);
        return new TanhSphericalDiffuse(eL, eR, w, c, integrator);
    }
    const std::string TANHSPHERICALDIFFUSE("TANHSPHERICALDIFFUSE");
    const bool registeredTanhSphericalDiffuse =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            TANHSPHERICALDIFFUSE, createTanhSphericalDiffuse);
}

#endif // TANHSPHERICALDIFFUSE_HPP
