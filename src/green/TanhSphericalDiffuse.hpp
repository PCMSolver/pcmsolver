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
#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>
#include "taylor.hpp"

// Boost general purpose includes
#include <boost/bind.hpp>
#include <boost/bind/placeholders.hpp>
#include <boost/function.hpp>
// Boost.Odeint includes
#include <boost/numeric/odeint/config.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>

#include "DerivativeTypes.hpp"
#include "DiagonalIntegratorFactory.hpp"
#include "DiagonalIntegrator.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"
#include "SphericalDiffuse.hpp"
#include "TanhDiffuse.hpp"
#include "Timer.hpp"
#include "LoggerInterface.hpp"

/*! \file TanhSphericalDiffuse.hpp
 *  \typedef TanhSphericalDiffuse
 *  \brief Green's function for a diffuse interface with spherical symmetry and tanh profile
 *  \author Hui Cao, Ville Weijo, Luca Frediani and Roberto Di Remigio
 *  \date 2010-2015
 */
typedef SphericalDiffuse<TanhDiffuse> TanhSphericalDiffuse;

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
        void l(size_t val) { l_ = val; }
        /*! Constructor from profile evaluator */
        LnTransformedRadial(const ProfileEvaluator & e) : eval_(e), l_(0) {}
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
 *  \struct Observer
 *  \brief reports progress of differential equation integrator
 *  \author Roberto Di Remigio
 *  \date 2015
 */
struct Observer
{
    void operator()(RadialFunction & f, const StateType & x, double r) const
    {
        /* Save grid points */
        f[0].push_back(r);
        /* Save function */
        f[1].push_back(x[0]);
        /* Save first derivative of function */
        f[2].push_back(x[1]);
    }
};

/*! \typedef IntegratorObserver
 *  \brief boosted-up function pointer to the observer operator()
 */
typedef boost::function<void (const StateType &, double)> IntegratorObserver;

template <>
inline void TanhSphericalDiffuse::computeZeta(const ProfileEvaluator & eval, const IntegratorParameters & params)
{
    namespace odeint = boost::numeric::odeint;
    odeint::bulirsch_stoer_dense_out<StateType> stepper_(params.eps_abs_, params.eps_rel_, params.factor_x_, params.factor_dxdt_);
    LnTransformedRadial system_(eval); /*! The system of first-order ODEs */
    StateType init_zeta_(2);           /*! Holds the initial conditions */
    Observer obs_;
    for (size_t L = 0; L < maxLGreen_; ++L) {
        // Create an empty RadialFunction
        RadialFunction f_;
        // Partial application of Observer::operator()
        IntegratorObserver observer_ = boost::bind(&Observer::operator(), obs_, f_, _1, _2);
        // Set angular momentum
        system_.l(L);
        // Set initial conditions
        init_zeta_[0] = L * std::log(params.r_0_);
        init_zeta_[1] = L / params.r_0_;
        odeint::integrate_adaptive(stepper_, system_, init_zeta_, params.r_0_, params.r_infinity_, params.observer_step_, observer_);
        zeta_.push_back(f_);
    }
}

template <>
inline void TanhSphericalDiffuse::computeOmega(const ProfileEvaluator & eval, const IntegratorParameters & params)
{
    namespace odeint = boost::numeric::odeint;
    odeint::bulirsch_stoer_dense_out<StateType> stepper_(params.eps_abs_, params.eps_rel_, params.factor_x_, params.factor_dxdt_);
    LnTransformedRadial system_(eval); /*! The system of first-order ODEs */
    StateType init_omega_(2);          /*! Holds the initial conditions */
    Observer obs_;
    for (size_t L = 0; L < maxLGreen_; ++L) {
        // Create an empty RadialFunction
        RadialFunction f_;
        // Partial application of Observer::operator()
        IntegratorObserver observer_ = boost::bind(&Observer::operator(), obs_, f_, _1, _2);
        // Set angular momentum
        system_.l(L);
        // Set initial conditions
        init_omega_[0] = -(L + 1) * std::log(params.r_infinity_);
        init_omega_[1] = -(L + 1) / params.r_infinity_;
        // Notice that we integrate BACKWARDS, so we pass -params.observer_step_ to integrate_adaptive
        odeint::integrate_adaptive(stepper_, system_, init_omega_, params.r_infinity_, params.r_0_, -params.observer_step_, observer_);
        omega_.push_back(f_);
    }
}

template <>
inline void TanhSphericalDiffuse::computeZetaAndOmega(const ProfileEvaluator & eval, const IntegratorParameters & params)
{
    namespace odeint = boost::numeric::odeint;
    odeint::bulirsch_stoer_dense_out<StateType> stepper_(params.eps_abs_, params.eps_rel_, params.factor_x_, params.factor_dxdt_);
    LnTransformedRadial system_(eval, maxLC_); /*! The system of first-order ODEs */
    Observer obs_;
    // Partial application of Observer::operator()
    IntegratorObserver observer_zeta_ = boost::bind(&Observer::operator(), obs_, zetaC_, _1, _2);
    // Set initial conditions for zeta
    StateType init_zeta_(2);
    init_zeta_[0] = maxLC_ * std::log(params.r_0_);
    init_zeta_[1] = maxLC_ / params.r_0_;
    // Integrate for zeta
    odeint::integrate_adaptive(stepper_, system_, init_zeta_, params.r_0_, params.r_infinity_, params.observer_step_, observer_zeta_);

    // Partial application of Observer::operator()
    IntegratorObserver observer_omega_ = boost::bind(&Observer::operator(), obs_, omegaC_, _1, _2);
    // Set initial conditions for omega
    StateType init_omega_(2);
    init_omega_[0] = -(maxLC_ + 1) * std::log(params.r_infinity_);
    init_omega_[1] = -(maxLC_ + 1) / params.r_infinity_;
    // Integrate for omega
    // Notice that we integrate BACKWARDS, so we pass -params.observer_step_ to integrate_adaptive
    odeint::integrate_adaptive(stepper_, system_, init_omega_, params.r_infinity_, params.r_0_, -params.observer_step_, observer_omega_);
}

template <>
inline void TanhSphericalDiffuse::initSphericalDiffuse()
{
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
    // Bind the TanhDiffuse::operator() to an evaluator function
    ProfileEvaluator eval_ = boost::bind(&TanhDiffuse::operator(), this->profile_, _1, _2, _3);

    LOG("Computing coefficient for the separation of the Coulomb singularity");
    timerON("TanhSphericalDiffuse::computeZetaAndOmega");
    computeZetaAndOmega(eval_, params_);
    timerOFF("TanhSphericalDiffuse::computeZetaAndOmega");

    LOG("Computing first radial solution");
    timerON("TanhSphericalDiffuse::computeZeta");
    computeZeta(eval_, params_);
    timerOFF("TanhSphericalDiffuse::computeZeta");

    LOG("Computing second radial solution");
    timerON("TanhSphericalDiffuse::computeOmega");
    computeOmega(eval_, params_);
    timerOFF("TanhSphericalDiffuse::computeOmega");
}


template <>
inline Numerical GreensFunction<Numerical, TanhDiffuse>::function(const Eigen::Vector3d & source,
                                        const Eigen::Vector3d & probe) const
{
    // To be implemented
    return 0.0;
}

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

/*! \file TanhSphericalDiffuse.hpp
 *  \class LnTransformedReducedRadial
 *  \brief system of ln-transformed first-order reduced radial differential equations
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  Provides a handle to the system of differential equations for the integrator.
 *  The dielectric profile comes in as a boost::function object.
 */
class LnTransformedReducedRadial
{
    private:
        /*! Dielectric profile function and derivative evaluation */
        ProfileEvaluator eval_;
        /*! Angular momentum */
        size_t l_;
    public:
        void l(size_t val) { l_ = val; }
        /*! Constructor from profile evaluator */
        LnTransformedReducedRadial(const ProfileEvaluator & e) : eval_(e) { l_ = 0; }
        /*! Provides a functor for the evaluation of the system
         *  of first-order ODEs needed by Boost.Odeint
         *  The second-order ODE and the system of first-order ODEs
         *  are reported in the manuscript.
         *  \param[in] xi state vector holding the function and its first derivative
         *  \param[out] dxidr state vector holding the first and second derivative
         *  \param[in] r position on the integration grid
         */
        void operator()(const StateType & xi, StateType & dxidr, const double r)
        {
            // Evaluate the dielectric profile
            double eps = 0.0, epsPrime = 0.0;
            eval_(eps, epsPrime, r);
            double gamma_epsilon = epsPrime / eps;
            // System of equations is defined here
            dxidr[0] = xi[1];
            dxidr[1] = -xi[1] * (xi[1] + gamma_epsilon * (1 - 1/r)) + l_ * (l_ + 1) / std::pow(r, 2);
        }
};

#endif // TANHSPHERICALDIFFUSE_HPP
