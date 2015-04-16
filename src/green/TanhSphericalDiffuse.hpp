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

#include <boost/bind.hpp>
#include <boost/bind/placeholders.hpp>
#include <boost/function.hpp>
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

/*! \file TanhSphericalDiffuse.hpp
 *  \typedef TanhSphericalDiffuse
 *  \brief Green's function for a diffuse interface with spherical symmetry and tanh profile
 *  \author Hui Cao, Ville Weijo, Luca Frediani and Roberto Di Remigio
 *  \date 2010-2015
 */
typedef SphericalDiffuse<TanhDiffuse> TanhSphericalDiffuse;

/*! \typedef StateType
 *  \brief state vector for the differential equation integrator
 */
typedef std::vector<double> StateType;
/*! \typedef ProfileEvaluator
 *  \brief boosted-up function pointer to the dielectric profile evaluation function
 */
typedef boost::function<void (double, double, double)> ProfileEvaluator;

/*! \file TanhSphericalDiffuse.hpp
 *  \class LogTransformedRadial
 *  \brief system of log-transformed first-order radial differential equations
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  Provides a handle to the system of differential equations for the integrator.
 *  The dielectric profile comes in as a boost::function object.
 */
class LogTransformedRadial
{
    private:
        /*! Angular momentum */
        size_t l_;
        /*! Dielectric profile function and derivative evaluation */
        ProfileEvaluator eval_;
    public:
        LogTransformedRadial( int l, const ProfileEvaluator & e)
            : l_(l), eval_(e) { }
        void operator()( const StateType & R , StateType & dRdr , const double r)
        {
            // Evaluate the dielectric profile
            double eps = 0.0, epsPrime = 0.0;
            eval_(eps, epsPrime, r);
            // System of equations is defined here
            dRdr[0] = R[1];
            dRdr[1] = -R[1]*(2.0/r + epsPrime/eps + R[1]) + l_ * (l_ + 1) / std::pow(r, 2);
        }
};

/*! \file TanhSphericalDiffuse.hpp
 *  \struct IntegratorObserver
 *  \brief reports progress of differential equation integrator
 *  \author Roberto Di Remigio
 *  \date 2015
 */
struct IntegratorObserver
{
    void operator()(const StateType & x, const double r) const
    {
        std::cout << r << "     " << x[0] << "      " << x[1] << std::endl;
    }
};

template <>
inline void TanhSphericalDiffuse::initSphericalDiffuse()
{
    namespace odeint = boost::numeric::odeint;
    // Parameters for the numerical solution of the radial differential equation
    /*! Maximum angular momentum in the final summation over Legendre polynomials to obtain G */
    size_t maxLGreen_ = 30;
    /*! Maximum angular momentum to obtain C(r, r'), needed to separate the Coulomb singularity */
    size_t maxLC_     = maxLGreen_ + 30;
    /*! Number of integration steps (determined by the adaptive integrator) */
    size_t nSteps_ = 0;
    // Bind the TanhDiffuse::operator() to an evaluator function
    ProfileEvaluator eval = boost::bind(&TanhDiffuse::operator(), this->profile_, _1, _2, _3);
    // Initialize the system of differential equations
    // Integrator parameters
    double eps_abs_     = 1.0e-8; /*! Absolute tolerance level */
    double eps_rel_     = 0.0;    /*! Relative tolerance level */
    double factor_x_    = 0.0;    /*! Weight of the state      */
    double factor_dxdt_ = 0.0;    /*! Weight of the state derivative */
    odeint::bulirsch_stoer_dense_out<StateType> stepper( eps_abs_, eps_rel_, factor_x_, factor_dxdt_ );
    for (size_t l = 0; l < maxLGreen_; ++l) {
        LogTransformedRadial system(l, eval);
    }
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

#endif // TANHSPHERICALDIFFUSE_HPP
