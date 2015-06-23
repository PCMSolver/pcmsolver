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

#ifndef INTERFACESIMPL_HPP
#define INTERFACESIMPL_HPP

#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <tuple>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

// Boost.Odeint includes
#include <boost/numeric/odeint.hpp>

#include "MathUtils.hpp"

namespace interfaces {
/*! \typedef StateType
 *  \brief state vector for the differential equation integrator
 */
typedef std::vector<double> StateType;

/*! \typedef RadialFunction
 *  \brief holds a solution to the radial equation: grid, function and first derivative
 */
typedef std::array<StateType, 3> RadialFunction;

/*! \typedef ProfileEvaluator
 *  \brief sort of a function pointer to the dielectric profile evaluation function
 */
typedef std::function< std::tuple<double, double>(const double) > ProfileEvaluator;

/*! \struct IntegratorParameters
 *  \brief holds parameters for the integrator
 */
struct IntegratorParameters
{
    /*! Absolute tolerance level */
    double eps_abs_     ;
    /*! Relative tolerance level */
    double eps_rel_     ;
    /*! Weight of the state      */
    double factor_x_    ;
    /*! Weight of the state derivative */
    double factor_dxdt_ ;
    /*! Lower bound of the integration interval */
    double r_0_         ;
    /*! Upper bound of the integration interval */
    double r_infinity_  ;
    /*! Time step between observer calls */
    double observer_step_;
    IntegratorParameters(double e_abs, double e_rel, double f_x, double f_dxdt, double r0, double rinf, double step)
    	: eps_abs_(e_abs), eps_rel_(e_rel), factor_x_(f_x),
    	factor_dxdt_(f_dxdt), r_0_(r0), r_infinity_(rinf), observer_step_(step) {}
};

/*! \class LnTransformedRadial
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
            std::tie(eps, epsPrime) = eval_(r);
            if (numericalZero(eps)) throw std::domain_error("Division by zero!");
            double gamma_epsilon = epsPrime / eps;
            // System of equations is defined here
            drhodr[0] = rho[1];
            drhodr[1] = -rho[1] * (rho[1] + 2.0/r + gamma_epsilon) + l_ * (l_ + 1) / std::pow(r, 2);
        }
};

/*! \brief reports progress of differential equation integrator
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
 *  \param[in, out] f RadialFunction whose contents have to be reversed
 *  \author Roberto Di Remigio
 *  \date 2015
 */
inline void reverse(RadialFunction & f)
{
    std::reverse(f[0].begin(), f[0].end());
    std::reverse(f[1].begin(), f[1].end());
    std::reverse(f[2].begin(), f[2].end());
}

/*! \brief reverse contents of a RadialFunction
 *  \param[in] f RadialFunction whose contents have to be printed
 *  \author Roberto Di Remigio
 *  \date 2015
 */
inline void writeRadialFunction(const RadialFunction & f, const std::string & fname)
{
    std::ofstream fout;
    fout.open(fname.c_str());
    fout << "#   r        f        df    "  << std::endl;
    int size = f[0].size();
    for (int i = 0; i < size; ++i) {
        fout << f[0][i] << "    " << f[1][i] << "      " << f[2][i] << std::endl;
    }
    fout.close();
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
    odeint::bulirsch_stoer_dense_out<StateType> stepper_(params.eps_abs_, params.eps_rel_, params.factor_x_, params.factor_dxdt_);
    // The system of first-order ODEs
    LnTransformedRadial system_(eval, L);
    // Holds the initial conditions
    StateType init_zeta_(2);
    // Set initial conditions
    init_zeta_[0] = L * std::log(params.r_0_);
    init_zeta_[1] = L / params.r_0_;
    odeint::integrate_adaptive(stepper_, system_, init_zeta_,
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
    odeint::bulirsch_stoer_dense_out<StateType> stepper_(params.eps_abs_, params.eps_rel_, params.factor_x_, params.factor_dxdt_);
    // The system of first-order ODEs
    LnTransformedRadial system_(eval, L);
    // Holds the initial conditions
    StateType init_omega_(2);
    // Set initial conditions
    init_omega_[0] = -(L + 1) * std::log(params.r_infinity_);
    init_omega_[1] = -(L + 1) / params.r_infinity_;
    // Notice that we integrate BACKWARDS, so we pass -params.observer_step_ to integrate_adaptive
    odeint::integrate_adaptive(stepper_, system_, init_omega_,
            params.r_infinity_, params.r_0_, -params.observer_step_,
            std::bind(observer, std::ref(f), _1, _2));
    // Reverse order of StateType-s in RadialFunction
    // this ensures that they are in ascending order (as later expected by linearInterpolation)
    reverse(f);
}
} // namespace interfaces

#endif // INTERFACESIMPL_HPP
