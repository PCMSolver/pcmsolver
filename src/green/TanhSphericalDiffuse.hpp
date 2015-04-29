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
#include <boost/numeric/odeint.hpp>
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
#include "GeneralPurpose.hpp"

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
        void operator()(const StateType & rho, StateType & drhodr, const double r)
        {
            // Evaluate the dielectric profile
            double eps = 0.0, epsPrime = 0.0;
            eval_(eps, epsPrime, r);
            if (numericalZero(eps)) throw std::invalid_argument("Division by zero!");
            double gamma_epsilon = epsPrime / eps;
            // System of equations is defined here
            drhodr[0] = rho[1];
            drhodr[1] = -rho[1] * (rho[1] + 2.0/r + gamma_epsilon) + l_ * (l_ + 1) / std::pow(r, 2);
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

/*! \brief reverse contents of a RadialFunction
 *  \param[in] f RadialFunction whose contents have to be reversed
 *  \author Roberto Di Remigio
 *  \date 2015
 */
inline void reverse(RadialFunction & f)
{
    std::reverse(f[0].begin(), f[0].end());
    std::reverse(f[1].begin(), f[1].end());
    std::reverse(f[2].begin(), f[2].end());
}

/*! \brief Calculates 1st radial solution, i.e. the one with r^l behavior
 *  \param[in]  L      angular momentum of the required solution
 *  \param[out] f      solution to the radial equation
 *  \param[in]  eval   dielectric profile evaluator function object
 *  \param[in]  params parameters for the integrator
 */
inline void computeZeta(int L, RadialFunction & f, const ProfileEvaluator & eval, const IntegratorParameters & params)
{
    using namespace std::placeholders;
    namespace odeint = boost::numeric::odeint;
    odeint::runge_kutta_fehlberg78<StateType> stepper_;
    // The system of first-order ODEs
    LnTransformedRadial system_(eval, L);
    // Holds the initial conditions
    StateType init_zeta_(2);
    // Set initial conditions
    init_zeta_[0] = L * std::log(params.r_0_);
    init_zeta_[1] = L / params.r_0_;
    odeint::integrate_const(stepper_, system_, init_zeta_,
            params.r_0_, params.r_infinity_, params.observer_step_,
            std::bind(observer, std::ref(f), _1, _2));
}

/*! \brief Calculates 2nd radial solution, i.e. the one with r^(-l-1) behavior
 *  \param[in]  L      angular momentum of the required solution
 *  \param[out] f      solution to the radial equation
 *  \param[in]  eval   dielectric profile evaluator function object
 *  \param[in]  params parameters for the integrator
 */
inline void computeOmega(int L, RadialFunction & f, const ProfileEvaluator & eval, const IntegratorParameters & params)
{
    using namespace std::placeholders;
    namespace odeint = boost::numeric::odeint;
    odeint::runge_kutta_fehlberg78<StateType> stepper_;
    // The system of first-order ODEs
    LnTransformedRadial system_(eval, L);
    // Holds the initial conditions
    StateType init_omega_(2);
    // Set initial conditions
    init_omega_[0] = -(L + 1) * std::log(params.r_infinity_);
    init_omega_[1] = -(L + 1) / params.r_infinity_;
    // Notice that we integrate BACKWARDS, so we pass -params.observer_step_ to integrate_adaptive
    odeint::integrate_const(stepper_, system_, init_omega_,
            params.r_infinity_, params.r_0_, -params.observer_step_,
            std::bind(observer, std::ref(f), _1, _2));
    // Reverse order of StateType-s in RadialFunction
    // this ensures that they are in ascending order (as later expected by linearInterpolation)
    reverse(f);
}

template <>
inline void TanhSphericalDiffuse::initSphericalDiffuse()
{
    using namespace std::placeholders;

    LOG("TanhSphericalDiffuse::initSphericalDiffuse");
    // Parameters for the numerical solution of the radial differential equation
    double eps_abs_     = 1.0e-10; /*! Absolute tolerance level */
    double eps_rel_     = 1.0e-06; /*! Relative tolerance level */
    double factor_x_    = 1.0;     /*! Weight of the state      */
    double factor_dxdt_ = 1.0;     /*! Weight of the state derivative */
    double r_0_         = 0.5;     /*! Lower bound of the integration interval */
    double r_infinity_  = profile_.center() + 100.0; /*! Upper bound of the integration interval */
    double observer_step_ = 1.0e-03; /*! Time step between observer calls */
    IntegratorParameters params_(eps_abs_, eps_rel_, factor_x_, factor_dxdt_, r_0_, r_infinity_, observer_step_);
    ProfileEvaluator eval_ = std::bind(&TanhDiffuse::operator(), this->profile_, _1, _2, _3);

    LOG("Computing coefficient for the separation of the Coulomb singularity");
    LOG("Computing first radial solution L = " + std::to_string(maxLC_));
    timerON("TanhSphericalDiffuse::computeZeta for coefficient");
    computeZeta(maxLC_, zetaC_, eval_, params_);
    timerOFF("TanhSphericalDiffuse::computeZeta for coefficient");
    LOG("DONE: Computing first radial solution L = " + std::to_string(maxLC_));

    LOG("Computing second radial solution L = " + std::to_string(maxLC_));
    timerON("TanhSphericalDiffuse::computeOmega");
    computeOmega(maxLC_, omegaC_, eval_, params_);
    timerOFF("TanhSphericalDiffuse::computeOmega");
    LOG("Computing second radial solution L = " + std::to_string(maxLC_));
    LOG("DONE: Computing coefficient for the separation of the Coulomb singularity");

    LOG("Computing radial solutions for Green's function");
    timerON("   Looping over angular momentum");
    for (int L = 0; L <= maxLGreen_; ++L) {
        // First radial solution
        LOG("Computing first radial solution L = " + std::to_string(L));
        timerON("TanhSphericalDiffuse::computeZeta");
        // Create an empty RadialFunction
        RadialFunction tmp_zeta_;
        computeZeta(L, tmp_zeta_, eval_, params_);
        zeta_.push_back(tmp_zeta_);
        timerOFF("TanhSphericalDiffuse::computeZeta");
        LOG("DONE: Computing first radial solution L = " + std::to_string(L));

        // Second radial solution
        LOG("Computing second radial solution L = " + std::to_string(L));
        timerON("TanhSphericalDiffuse::computeOmega");
        // Create an empty RadialFunction
        RadialFunction tmp_omega_;
        computeOmega(L, tmp_omega_, eval_, params_);
        omega_.push_back(tmp_omega_);
        timerOFF("TanhSphericalDiffuse::computeOmega");
        LOG("DONE: Computing second radial solution L = " + std::to_string(L));
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
inline double TanhSphericalDiffuse::functionSummation(int L, double r1, double r2, double cos_gamma, double Cr12) const
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

template <>
inline double TanhSphericalDiffuse::Coulomb(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const
{
    double r1  = source.norm();
    double r2  = probe.norm();
    double r12 = (source - probe).norm();
    double cos_gamma = source.dot(probe) / (r1 * r2);

    // Obtain coefficient for the separation of the Coulomb singularity
    return (1.0 / (this->coefficient(r1, r2) * r12));
}

template <>
inline double TanhSphericalDiffuse::imagePotential(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const
{
    double r1  = source.norm();
    double r2  = probe.norm();
    double r12 = (source - probe).norm();
    double cos_gamma = source.dot(probe) / (r1 * r2);

    // Obtain coefficient for the separation of the Coulomb singularity
    double Cr12 = this->coefficient(r1, r2);

    double gr12 = 0.0;
    for (int L = 0; L <= maxLGreen_; ++L) {
        gr12 += this->functionSummation(L, r1, r2, cos_gamma, Cr12);
    }

    return gr12;
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
inline Eigen::Vector3d TanhSphericalDiffuse::coefficientGradient(const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d coeff_grad = Eigen::Vector3d::Zero();

    double r1  = p1.norm();
    double r2  = p2.norm();
    double r2_2 = std::pow(r2, 2);
    double r2_3 = std::pow(r2, 3);
    double r12 = (p1 - p2).norm();
    double cos_gamma = p1.dot(p2) / (r1 * r2);
    double cos_gamma_2 = std::pow(cos_gamma, 2);

    double pl_x = boost::math::legendre_p(maxLC_, cos_gamma);
    double pl_first_x = -boost::math::legendre_p(maxLC_, 1, cos_gamma)
                                 / (std::sqrt(1.0 - std::pow(cos_gamma, 2)));

    double eps_r2 = 0.0, epsPrime_r2 = 0.0;
    this->profile_(eps_r2, epsPrime_r2, r2);

    /* Value of zetaC_ at point with index 1 */
    double zeta1  = linearInterpolation(r1, zetaC_[0], zetaC_[1]);
    /* Value of zetaC_ at point with index 2 */
    double zeta2  = linearInterpolation(r2, zetaC_[0], zetaC_[1]);
    /* Value of omegaC_ at point with index 1 */
    double omega1 = linearInterpolation(r1, omegaC_[0], omegaC_[1]);
    /* Value of omegaC_ at point with index 2 */
    double omega2 = linearInterpolation(r2, omegaC_[0], omegaC_[1]);

    /* Value of derivative of zetaC_ at point with index 2 */
    double d_zeta2  = linearInterpolation(r2, zetaC_[0], zetaC_[2]);
    /* Value of derivative of omegaC_ at point with index 2 */
    double d_omega2 = linearInterpolation(r2, omegaC_[0], omegaC_[2]);

    /*! Value of second derivative of zetaC_ at point with index 2 */
    double d2_zeta2 = maxLC_ * (maxLC_ + 1) / r2_2 - std::pow(d_zeta2, 2) - (2.0 / r2 + epsPrime_r2 / eps_r2) * d_zeta2;
    /*! Value of second derivative of omegaC_ at point with index 2 */
    double d2_omega2 = maxLC_ * (maxLC_ + 1) / r2_2 - std::pow(d_omega2, 2) - (2.0 / r2 + epsPrime_r2 / eps_r2) * d_omega2;

    double a = (epsPrime_r2 * r2 + 2.0 * eps_r2) * pl_x / (eps_r2 * r2);
    double b = (d_zeta2 * (d_zeta2 - d_omega2) + d2_zeta2 - d2_omega2) * pl_x / (d_zeta2 - d_omega2);
    double c = (d_omega2 * (d_zeta2 - d_omega2) + d2_zeta2 - d2_omega2) * pl_x / (d_zeta2 - d_omega2);

    double x1 = p1(0);
    double y1 = p1(1);
    double z1 = p1(2);
    double x2 = p2(0);
    double y2 = p2(1);
    double z2 = p2(2);

    Eigen::Vector3d tmp_grad = Eigen::Vector3d::Zero();
    if (r1 < r2) {
        double expFact = std::exp(zeta1 -zeta2) * (2*maxLC_ + 1);
        tmp_grad(0) = -expFact*x2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + b - (x1*r2*r2/x2-r12)*pl_first_x/
                      (r1*r2*r2));
        tmp_grad(1) = -expFact*y2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + b - (y1*r2*r2/y2-r12)*pl_first_x/
                      (r1*r2*r2));
        tmp_grad(2) = -expFact*z2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + b - (z1*r2*r2/z2-r12)*pl_first_x/
                      (r1*r2*r2));

        double powl = std::pow(r1/r2, maxLC_);
        coeff_grad(0) = powl*(pl_x*p2(0)*(-maxLC_-1)+pl_first_x*(p1(0)*r2_2 -p2(0)*(p1.dot(
                                    p2)))/(r1*r2))/(r2_3*tmp_grad(0));
        coeff_grad(1) = powl*(pl_x*p2(1)*(-maxLC_-1)+pl_first_x*(p1(1)*r2_2-p2(1)*(p1.dot(
                                    p2)))/(r1*r2))/(r2_3*tmp_grad(1));
        coeff_grad(2) = powl*(pl_x*p2(2)*(-maxLC_-1)+pl_first_x*(p1(2)*r2_2-p2(2)*(p1.dot(
                                    p2)))/(r1*r2))/(r2_3*tmp_grad(2));
    } else {
        double expFact = exp(omega1 - omega2) * (2*maxLC_ + 1);
        tmp_grad(0) = -expFact*x2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + c - (x1*r2*r2/x2-r12)*pl_first_x/
                      (r1*r2*r2));
        tmp_grad(1) = -expFact*y2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + c - (y1*r2*r2/y2-r12)*pl_first_x/
                      (r1*r2*r2));
        tmp_grad(2) = -expFact*z2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + c - (z1*r2*r2/z2-r12)*pl_first_x/
                      (r1*r2*r2));

        double powl = std::pow(r2/r1, maxLC_);
        coeff_grad(0) = powl*(pl_x*p2(0)*maxLC_+pl_first_x*(p1(0)*r2_2 -p2(0)*(p1.dot(p2)))/
                                (r1*r2))/(r1*r2_2*tmp_grad(0));
        coeff_grad(1) = powl*(pl_x*p2(1)*maxLC_+pl_first_x*(p1(1)*r2_2 -p2(1)*(p1.dot(p2)))/
                                (r1*r2))/(r1*r2_2*tmp_grad(1));
        coeff_grad(2) = powl*(pl_x*p2(2)*maxLC_+pl_first_x*(p1(2)*r2_2 -p2(2)*(p1.dot(p2)))/
                                (r1*r2))/(r1*r2_2*tmp_grad(2));
    }

    return coeff_grad;
}

template <>
inline Eigen::Vector3d TanhSphericalDiffuse::functionSummationGradient(int L, const Eigen::Vector3d & p1,
        const Eigen::Vector3d & p2, double Cr12) const
{
    Eigen::Vector3d gr12_grad = Eigen::Vector3d::Zero();
    return gr12_grad;
}

template <>
inline Numerical GreensFunction<Numerical, TanhDiffuse>::derivativeProbe(const Eigen::Vector3d & normal_p2,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    double eps_r2 = 0.0, epsPrime_r2 = 0.0;
    this->profile_(eps_r2, epsPrime_r2, p2.norm());
    Eigen::Vector3d grad = Eigen::Vector3d::Zero();

    // To be implemented
    return (eps_r2 * grad.dot(normal_p2));
}

/*
Eigen::Vector3d SphericalInterface::converged_deri_gf(const Eigen::VectorXd & p1,
        const Eigen::VectorXd & p2, double *plx,
        double *dplx) const
{
    Eigen::Vector3d greenGradient(0.0, 0.0, 0.0);
    Eigen::Vector3d Cr12Gradient(0.0, 0.0, 0.0);
    double r1 = p1.norm();
    double r2 = p2.norm();
    double eps2, deps2;
    profile(&eps2, &deps2, r2);

    double pl = plx[maxLEpsilon_];
    double dpl = dplx[maxLEpsilon_];

    Cr12Gradient.Zero();
    greenGradient = greenfunc_der(p1, p2, Cr12Gradient, radialC1_, radialC2_,
                                  pl, dpl, maxLEpsilon_, 0);
    if (r1 < r2) {
        double powl = std::pow(r1/r2, maxLEpsilon_);
        Cr12Gradient(0) = powl*(pl*p2(0)*(-maxLEpsilon_-1)+dpl*(p1(0)*r2*r2-p2(0)*(p1.dot(
                                    p2)))/(r1*r2))/(r2*r2*r2*greenGradient(0));
        Cr12Gradient(1) = powl*(pl*p2(1)*(-maxLEpsilon_-1)+dpl*(p1(1)*r2*r2-p2(1)*(p1.dot(
                                    p2)))/(r1*r2))/(r2*r2*r2*greenGradient(1));
        Cr12Gradient(2) = powl*(pl*p2(2)*(-maxLEpsilon_-1)+dpl*(p1(2)*r2*r2-p2(2)*(p1.dot(
                                    p2)))/(r1*r2))/(r2*r2*r2*greenGradient(2));
    } else {
        double powl = std::pow(r2/r1, maxLEpsilon_);
        Cr12Gradient(0) = powl*(pl*p2(0)*maxLEpsilon_+dpl*(p1(0)*r2*r2-p2(0)*(p1.dot(p2)))/
                                (r1*r2))/(r1*r2*r2*greenGradient(0));
        Cr12Gradient(1) = powl*(pl*p2(1)*maxLEpsilon_+dpl*(p1(1)*r2*r2-p2(1)*(p1.dot(p2)))/
                                (r1*r2))/(r1*r2*r2*greenGradient(1));
        Cr12Gradient(2) = powl*(pl*p2(2)*maxLEpsilon_+dpl*(p1(2)*r2*r2-p2(2)*(p1.dot(p2)))/
                                (r1*r2))/(r1*r2*r2*greenGradient(2));
    }

    greenGradient.Zero();

    for (int l = 0; l < maxLGreen_; ++l) {
        Eigen::Vector3d gl = greenfunc_der(p1, p2, Cr12Gradient, radialG1_[l],
                                           radialG2_[l], plx[l], dplx[l], l, 1);
        greenGradient += gl;
    }

    double r12 = (p1-p2).norm();
    Eigen::Vector3d gr_d;
    gr_d.array() = (p1.array() - p2.array()) /
                   (Cr12Gradient.array() * r12 * r12 * r12) + greenGradient.array();
    gr_d *= -1;
    return gr_d;
}

// greenGradient = greenfunc_der(p1, p2, Cr12Gradient, radialC1_, radialC2_, pl, dpl, maxLEpsilon_, 0);

Eigen::Vector3d SphericalInterface::greenfunc_der(const Eigen::Vector3d & p1,
        const Eigen::Vector3d & p2, const Eigen::Vector3d & Cr12,
        const FuncGrid & f1, const FuncGrid & f2, double plx,
        double dplx, int l, int flagCr12) const
{
    double eps2, deps2;
    double prod = p1.dot(p2);
    double r1 = p1.norm();
    double r2 = p2.norm();
    double r2_3 = r2 * r2 * r2;
    profile(&eps2, &deps2, r2);

    //calculate Legendre function of L and x--cos_angle

    int idx1 = int((r1-rBegin_) / hStep_) - 1;
    int idx2 = int((r2-rBegin_) / hStep_) - 1;

    double delta1 = r1 - grid_(idx1);
    double delta2 = r2 - grid_(idx2);


    double u11 = f1(0, idx1) + (f1(0, idx1 + 1) - f1(0, idx1)) * delta1 / hStep_;
    double u12 = f1(0, idx2) + (f1(0, idx2 + 1) - f1(0, idx2)) * delta2 / hStep_;
    double u21 = f2(0, idx1) + (f2(0, idx1 + 1) - f2(0, idx1)) * delta1 / hStep_;
    double u22 = f2(0, idx2) + (f2(0, idx2 + 1) - f2(0, idx2)) * delta2 / hStep_;

    double d11 = f1(1, idx1) + (f1(1, idx1 + 1) - f1(1, idx1)) * delta1 / hStep_;
    double d12 = f1(1, idx2) + (f1(1, idx2 + 1) - f1(1, idx2)) * delta2 / hStep_;
    double d21 = f2(1, idx1) + (f2(1, idx1 + 1) - f2(1, idx1)) * delta1 / hStep_;
    double d22 = f2(1, idx2) + (f2(1, idx2 + 1) - f2(1, idx2)) * delta2 / hStep_;

    double d2u12 = l * (l + 1) / (r2*r2) - d12 * d12 - (2.0 / r2 + deps2 / eps2) * d12;
    double d2u22 = l * (l + 1) / (r2*r2) - d22 * d22 - (2.0 / r2 + deps2 / eps2) * d22;
    Eigen::Vector3d g12_deri;

    double x1 = p1(0);
    double y1 = p1(1);
    double z1 = p1(2);
    double x2 = p2(0);
    double y2 = p2(1);
    double z2 = p2(2);

    double a = (deps2*r2+2.0*eps2)*plx/(eps2*r2);
    double b = (d12*(d12-d22)+d2u12-d2u22)*plx/(d12-d22);
    double c = (d22*(d12-d22)+d2u12-d2u22)*plx/(d12-d22);
    if (r1<r2) {
        double expFact = exp(u11-u12)*(2*l+1);
        g12_deri(0) = -expFact*x2/((d12-d22)*r2_3*eps2) * (a + b - (x1*r2*r2/x2-prod)*dplx/
                      (r1*r2*r2));
        g12_deri(1) = -expFact*y2/((d12-d22)*r2_3*eps2) * (a + b - (y1*r2*r2/y2-prod)*dplx/
                      (r1*r2*r2));
        g12_deri(2) = -expFact*z2/((d12-d22)*r2_3*eps2) * (a + b - (z1*r2*r2/z2-prod)*dplx/
                      (r1*r2*r2));
        if (flagCr12 == 1) {
            g12_deri(0) -= (std::pow(r1/r2,
                                l)*(plx*x2*(-l-1)+dplx*(x1*r2*r2-x2*prod)/(r1*r2))/(r2_3*Cr12(0)));
            g12_deri(1) -= (std::pow(r1/r2,
                                l)*(plx*y2*(-l-1)+dplx*(y1*r2*r2-y2*prod)/(r1*r2))/(r2_3*Cr12(1)));
            g12_deri(2) -= (std::pow(r1/r2,
                                l)*(plx*z2*(-l-1)+dplx*(z1*r2*r2-z2*prod)/(r1*r2))/(r2_3*Cr12(2)));
        }
    } else {
        double expFact = exp(u21-u22)*(2*l+1);
        g12_deri(0) = -expFact*x2/((d12-d22)*r2_3*eps2) * (a + c - (x1*r2*r2/x2-prod)*dplx/
                      (r1*r2*r2));
        g12_deri(1) = -expFact*y2/((d12-d22)*r2_3*eps2) * (a + c - (y1*r2*r2/y2-prod)*dplx/
                      (r1*r2*r2));
        g12_deri(2) = -expFact*z2/((d12-d22)*r2_3*eps2) * (a + c - (z1*r2*r2/z2-prod)*dplx/
                      (r1*r2*r2));
        if (flagCr12 == 1) {
            g12_deri(0) -= (std::pow(r2/r1,
                                l)*(plx*x2*l+dplx*(x1*r2*r2-x2*prod)/(r1*r2))/(r1*r2*r2*Cr12(0)));
            g12_deri(1) -= (std::pow(r2/r1,
                                l)*(plx*y2*l+dplx*(y1*r2*r2-y2*prod)/(r1*r2))/(r1*r2*r2*Cr12(1)));
            g12_deri(2) -= (std::pow(r2/r1,
                                l)*(plx*z2*l+dplx*(z1*r2*r2-z2*prod)/(r1*r2))/(r1*r2*r2*Cr12(2)));
        }
    }
    return g12_deri;
}
*/

template <>
inline Numerical GreensFunction<Numerical, TanhDiffuse>::derivativeSource(const Eigen::Vector3d & normal_p1,
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
