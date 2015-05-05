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
#include "InterfacesDef.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"
#include "SphericalDiffuse.hpp"
#include "TanhDiffuse.hpp"
#include "Timer.hpp"
#include "LoggerInterface.hpp"

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
    double r_0_         = 0.5;     /*! Lower bound of the integration interval */
    double r_infinity_  = profile_.center() + 200.0; /*! Upper bound of the integration interval */
    double observer_step_ = 1.0e-03; /*! Time step between observer calls */
    IntegratorParameters params_(r_0_, r_infinity_, observer_step_);
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

/*! Calculates the Green's function given a pair of points */
template <>
inline double GreensFunction<Numerical, TanhDiffuse>::function(const Eigen::Vector3d & source,
                                        const Eigen::Vector3d & probe) const
{
    Numerical sp[3], pp[3];
    sp[0] = source(0); sp[1] = source(1); sp[2] = source(2);
    pp[0] = probe(0);  pp[1] = probe(1);  pp[2] = probe(2);
    return this->operator()(sp, pp);
}

template <>
inline double TanhSphericalDiffuse::coefficientCoulomb(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const
{
    double r1  = source.norm();
    double r2  = probe.norm();

    // Obtain coefficient for the separation of the Coulomb singularity
    return this->coefficient(r1, r2);
}

template <>
inline double TanhSphericalDiffuse::Coulomb(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const
{
    double r1  = source.norm();
    double r2  = probe.norm();
    double r12 = (source - probe).norm();

    // Obtain coefficient for the separation of the Coulomb singularity
    return (1.0 / (this->coefficient(r1, r2) * r12));
}

template <>
inline double TanhSphericalDiffuse::imagePotential(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const
{
    double r1  = source.norm();
    double r2  = probe.norm();
    double cos_gamma = source.dot(probe) / (r1 * r2);

    // Obtain coefficient for the separation of the Coulomb singularity
    double Cr12 = this->coefficient(r1, r2);

    double gr12 = 0.0;
    for (int L = 0; L <= maxLGreen_; ++L) {
        gr12 += this->functionSummation(L, r1, r2, cos_gamma, Cr12);
    }

    return gr12;
}

template <>
inline double TanhSphericalDiffuse::CoulombDerivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    using namespace std::placeholders;
    return threePointStencil(std::bind(&TanhSphericalDiffuse::Coulomb, this, _1, _2),
            p2, p1, direction, this->delta_);
}

template <>
inline double TanhSphericalDiffuse::imagePotentialDerivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    using namespace std::placeholders;
    return threePointStencil(std::bind(&TanhSphericalDiffuse::imagePotential, this, _1, _2),
                p2, p1, direction, this->delta_);
}

template <>
inline double GreensFunction<Numerical, TanhDiffuse>::derivativeProbe(const Eigen::Vector3d & /* normal_p2 */,
            const Eigen::Vector3d & /* p1 */, const Eigen::Vector3d & /* p2 */) const
{
    /*
    Eigen::Vector3d deltaPlus  = p2 + normal_p2 * this->delta_ / normal_p2.norm();
    Eigen::Vector3d deltaMinus = p2 - normal_p2 * this->delta_ / normal_p2.norm();
    double CoulombPlus = this->Coulomb(p1, deltaPlus);
    double CoulombMinus = this->Coulomb(p1, deltaMinus);
    double CoulombDeriv = (CoulombPlus - CoulombMinus) / (2.0 * this->delta_);
    double imagePlus = this->imagePotential(p1, deltaPlus);
    double imageMinus = this->imagePotential(p1, deltaMinus);
    double imageDeriv = (imagePlus - imageMinus) / (2.0 * this->delta_);

    return (CoulombDeriv + imageDeriv);
    */
    return 0.0;
}

template <>
inline double GreensFunction<Numerical, TanhDiffuse>::derivativeSource(const Eigen::Vector3d & /* normal_p1 */,
            const Eigen::Vector3d & /* p1 */, const Eigen::Vector3d & /* p2 */) const
{
    /*
    Eigen::Vector3d deltaPlus  = p1 + normal_p1 * this->delta_ / normal_p1.norm();
    Eigen::Vector3d deltaMinus = p1 - normal_p1 * this->delta_ / normal_p1.norm();
    double CoulombPlus = this->Coulomb(deltaPlus, p2);
    double CoulombMinus = this->Coulomb(deltaMinus, p2);
    double CoulombDeriv = (CoulombPlus - CoulombMinus) / (2.0 * this->delta_);
    double imagePlus = this->imagePotential(deltaPlus, p2);
    double imageMinus = this->imagePotential(deltaMinus, p2);
    double imageDeriv = (imagePlus - imageMinus) / (2.0 * this->delta_);

    return (CoulombDeriv + imageDeriv);
    */
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

/*
        // This should be the semi-analytic implementation...
        // But I get NaN-s and inf-s out of it, so I switched to
        // the good old three-point stencil.
template <>
inline double TanhSphericalDiffuse::derivative(const Eigen::Vector3d & direction,
                 const Eigen::Vector3d & p1, const Eigen::Vector3d & p2)
{
        double eps_r2 = 0.0, epsPrime_r2 = 0.0;
        this->profile_(eps_r2, epsPrime_r2, p2.norm());

        Eigen::Vector3d Cr12_grad = this->coefficientGradient(p1, p2);
        Eigen::Vector3d grad = Eigen::Vector3d::Zero();
        for (int L = 0; L < maxLGreen_; ++L) {
            grad += this->functionSummationGradient(L, p1, p2, Cr12_grad);
        }

        double r12 = (p1 - p2).norm();
        Eigen::Vector3d gr_d = Eigen::Vector3d::Zero();
        gr_d.array() = (p1.array() - p2.array()) /
	                   (Cr12_grad.array() * std::pow(r12, 3)) + grad.array();
        gr_d *= -1;

        return (eps_r2 * grad.dot(direction));
}

template <>
inline Eigen::Vector3d TanhSphericalDiffuse::coefficientGradient(const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d coeff_grad = Eigen::Vector3d::Zero();

    double r1  = p1.norm();
    double r2  = p2.norm();
    double r2_2 = std::pow(r2, 2);
    double r2_3 = std::pow(r2, 3);
    double prod = p1.dot(p2);
    double cos_gamma = p1.dot(p2) / (r1 * r2);

    double pl_x = boost::math::legendre_p(maxLC_, cos_gamma);
    double pl_first_x = -boost::math::legendre_p(maxLC_, 1, cos_gamma)
                                 / (std::sqrt(1.0 - std::pow(cos_gamma, 2)));

    double eps_r2 = 0.0, epsPrime_r2 = 0.0;
    this->profile_(eps_r2, epsPrime_r2, r2);

    // Value of zetaC_ at point with index 1
    double zeta1  = linearInterpolation(r1, zetaC_[0], zetaC_[1]);
    // Value of zetaC_ at point with index 2/
    double zeta2  = linearInterpolation(r2, zetaC_[0], zetaC_[1]);
    // Value of omegaC_ at point with index 1
    double omega1 = linearInterpolation(r1, omegaC_[0], omegaC_[1]);
    // Value of omegaC_ at point with index 2
    double omega2 = linearInterpolation(r2, omegaC_[0], omegaC_[1]);

    // Value of derivative of zetaC_ at point with index 2
    double d_zeta2  = linearInterpolation(r2, zetaC_[0], zetaC_[2]);
    // Value of derivative of omegaC_ at point with index 2
    double d_omega2 = linearInterpolation(r2, omegaC_[0], omegaC_[2]);

    // Value of second derivative of zetaC_ at point with index 2
    double d2_zeta2 = maxLC_ * (maxLC_ + 1) / r2_2 - std::pow(d_zeta2, 2) - (2.0 / r2 + epsPrime_r2 / eps_r2) * d_zeta2;
    // Value of second derivative of omegaC_ at point with index 2
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

    double factor_x = (x1 * r2_2 - x2 * prod) / (r1*r2);
    double factor_y = (y1 * r2_2 - y2 * prod) / (r1*r2);
    double factor_z = (z1 * r2_2 - z2 * prod) / (r1*r2);

    Eigen::Vector3d tmp_grad = Eigen::Vector3d::Zero();
    if (r1 < r2) {
        double expFact = std::exp(zeta1 -zeta2) * (2*maxLC_ + 1);
        tmp_grad(0) = -expFact*x2/((d_zeta2-d_omega2)*r2_3*eps_r2) *
            (a + b - (x1*r2_2/(x2-prod)) * pl_first_x / (r1*r2_2));
        tmp_grad(1) = -expFact*y2/((d_zeta2-d_omega2)*r2_3*eps_r2) *
            (a + b - (y1*r2_2/(y2-prod))*pl_first_x / (r1*r2_2));
        tmp_grad(2) = -expFact*z2/((d_zeta2-d_omega2)*r2_3*eps_r2) *
            (a + b - (z1*r2_2/(z2-prod))*pl_first_x/(r1*r2_2));

        double powl = std::pow(r1 / r2, maxLC_);
        coeff_grad(0) = powl * (pl_x * x2 * (-maxLC_-1) + pl_first_x * factor_x) / (r2_3 * tmp_grad(0));
        coeff_grad(1) = powl * (pl_x * y2 * (-maxLC_-1) + pl_first_x * factor_y) / (r2_3 * tmp_grad(1));
        coeff_grad(2) = powl * (pl_x * z2 * (-maxLC_-1) + pl_first_x * factor_z) / (r2_3 * tmp_grad(2));
    } else {
        double expFact = exp(omega1 - omega2) * (2*maxLC_ + 1);
        tmp_grad(0) = -expFact*x2/((d_zeta2-d_omega2)*r2_3*eps_r2) *
            (a + c - (x1*r2_2/(x2 - prod)) * pl_first_x / (r1*r2_2));
        tmp_grad(1) = -expFact*y2/((d_zeta2-d_omega2)*r2_3*eps_r2) *
            (a + c - (y1*r2_2/(y2 - prod)) * pl_first_x / (r1*r2_2));
        tmp_grad(2) = -expFact*z2/((d_zeta2-d_omega2)*r2_3*eps_r2)
            * (a + c - (z1*r2_2/(z2 - prod)) * pl_first_x / (r1*r2_2));

        double powl = std::pow(r2 / r1, maxLC_);
        coeff_grad(0) = powl * (pl_x * x2 * maxLC_ + pl_first_x * factor_x) / (r1 * r2_2 * tmp_grad(0));
        coeff_grad(1) = powl * (pl_x * y2 * maxLC_ + pl_first_x * factor_y) / (r1 * r2_2 * tmp_grad(1));
        coeff_grad(2) = powl * (pl_x * z2 * maxLC_ + pl_first_x * factor_z) / (r1 * r2_2 * tmp_grad(2));
    }

    return coeff_grad;
}

template <>
inline Eigen::Vector3d TanhSphericalDiffuse::functionSummationGradient(int L, const Eigen::Vector3d & p1,
        const Eigen::Vector3d & p2, const Eigen::Vector3d & Cr12_grad) const
{
    Eigen::Vector3d gr12_grad = Eigen::Vector3d::Zero();

    double r1  = p1.norm();
    double r2  = p2.norm();
    double r2_2 = std::pow(r2, 2);
    double r2_3 = std::pow(r2, 3);
    double prod = p1.dot(p2);
    double cos_gamma = p1.dot(p2) / (r1 * r2);

    double pl_x = boost::math::legendre_p(maxLC_, cos_gamma);
    double pl_first_x = -boost::math::legendre_p(maxLC_, 1, cos_gamma)
                                 / (std::sqrt(1.0 - std::pow(cos_gamma, 2)));

    double eps_r2 = 0.0, epsPrime_r2 = 0.0;
    this->profile_(eps_r2, epsPrime_r2, r2);

    // Value of zetaC_ at point with index 1
    double zeta1  = linearInterpolation(r1, zetaC_[0], zetaC_[1]);
    // Value of zetaC_ at point with index 2
    double zeta2  = linearInterpolation(r2, zetaC_[0], zetaC_[1]);
    // Value of omegaC_ at point with index 1
    double omega1 = linearInterpolation(r1, omegaC_[0], omegaC_[1]);
    // Value of omegaC_ at point with index 2
    double omega2 = linearInterpolation(r2, omegaC_[0], omegaC_[1]);

    // Value of derivative of zetaC_ at point with index 2
    double d_zeta2  = linearInterpolation(r2, zetaC_[0], zetaC_[2]);
    // Value of derivative of omegaC_ at point with index 2
    double d_omega2 = linearInterpolation(r2, omegaC_[0], omegaC_[2]);

    // Value of second derivative of zetaC_ at point with index 2
    double d2_zeta2 = maxLC_ * (maxLC_ + 1) / r2_2 - std::pow(d_zeta2, 2) - (2.0 / r2 + epsPrime_r2 / eps_r2) * d_zeta2;
    // Value of second derivative of omegaC_ at point with index 2
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

	if (r1 < r2) {
		double expFact = exp(zeta1-zeta2)*(2*L+1);
		gr12_grad(0) = -expFact*x2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + b - (x1*r2*r2/(x2-prod))*pl_first_x/(r1*r2*r2));
		gr12_grad(1) = -expFact*y2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + b - (y1*r2*r2/(y2-prod))*pl_first_x/(r1*r2*r2));
		gr12_grad(2) = -expFact*z2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + b - (z1*r2*r2/(z2-prod))*pl_first_x/(r1*r2*r2));
		gr12_grad(0) -= (pow(r1/r2,L)*(pl_x*x2*(-L-1)+pl_first_x*(x1*r2*r2-x2*prod)/(r1*r2))/(r2_3*Cr12_grad(0)));
		gr12_grad(1) -= (pow(r1/r2,L)*(pl_x*y2*(-L-1)+pl_first_x*(y1*r2*r2-y2*prod)/(r1*r2))/(r2_3*Cr12_grad(1)));
		gr12_grad(2) -= (pow(r1/r2,L)*(pl_x*z2*(-L-1)+pl_first_x*(z1*r2*r2-z2*prod)/(r1*r2))/(r2_3*Cr12_grad(2)));
	} else {
		double expFact = exp(omega1-omega2)*(2*L+1);
		gr12_grad(0) = -expFact*x2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + c - (x1*r2*r2/(x2-prod))*pl_first_x/(r1*r2*r2));
		gr12_grad(1) = -expFact*y2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + c - (y1*r2*r2/(y2-prod))*pl_first_x/(r1*r2*r2));
		gr12_grad(2) = -expFact*z2/((d_zeta2-d_omega2)*r2_3*eps_r2) * (a + c - (z1*r2*r2/(z2-prod))*pl_first_x/(r1*r2*r2));
		gr12_grad(0) -= (pow(r2/r1,L)*(pl_x*x2*L+pl_first_x*(x1*r2*r2-x2*prod)/(r1*r2))/(r1*r2*r2*Cr12_grad(0)));
		gr12_grad(1) -= (pow(r2/r1,L)*(pl_x*y2*L+pl_first_x*(y1*r2*r2-y2*prod)/(r1*r2))/(r1*r2*r2*Cr12_grad(1)));
		gr12_grad(2) -= (pow(r2/r1,L)*(pl_x*z2*L+pl_first_x*(z1*r2*r2-z2*prod)/(r1*r2))/(r1*r2*r2*Cr12_grad(2)));
	}

    return gr12_grad;
}
*/

#endif // TANHSPHERICALDIFFUSE_HPP
