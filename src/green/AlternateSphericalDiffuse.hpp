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

#ifndef ALTERNATESPHERICALDIFFUSE_HPP
#define ALTERNATESPHERICALDIFFUSE_HPP

#include <array>
#include <cmath>
#include <functional>
#include <iosfwd>
#include <tuple>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

// Has to be included here
#include "InterfacesImpl.hpp"
// Boost.Math includes
#include <boost/math/special_functions/legendre.hpp>

#include "GreensFunction.hpp"
#include "LoggerInterface.hpp"
#include "MathUtils.hpp"
#include "Timer.hpp"

/*! \file AlternateSphericalDiffuse.hpp
 *  \class AlternateSphericalDiffuse
 *  \brief Green's function for a diffuse interface with spherical symmetry, with alternative separation of Coulomb singularity
 *  \author Hui Cao, Ville Weijo, Luca Frediani and Roberto Di Remigio
 *  \date 2010-2015
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation of S and D
 *  \tparam ProfilePolicy functional form of the diffuse layer
 *
 *  This class is general, in the sense that no specific dielectric profile has been set in its definition.
 *  In principle any profile that can be described by:
 *  1. a left-side dielectric constant;
 *  2. a right-side dielectric constant;
 *  3. an interface layer width;
 *  4. an interface layer center
 *  can be used to define a new diffuse interface with spherical symmetry.
 *  The origin of the dielectric sphere can be changed by means of the constructor.
 *  The solution of the differential equation defining the Green's function is **always**
 *  performed assuming that the dielectric sphere is centered in the origin of the coordinate
 *  system. Whenever the public methods are invoked to "sample" the Green's function
 *  at a pair of points, a translation of the sampling points is performed first.
 */

template <typename IntegratorPolicy,
          typename ProfilePolicy>
class AlternateSphericalDiffuse : public GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy,
                                               AlternateSphericalDiffuse<IntegratorPolicy, ProfilePolicy> >
{
public:
    /*! Constructor for a one-layer interface
     * \param[in] e1 left-side dielectric constant
     * \param[in] e2 right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] c center of the diffuse layer
     * \param[in] o center of the sphere
     */
    AlternateSphericalDiffuse(double e1, double e2, double w, double c, const Eigen::Vector3d & o, int l)
        : GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy, AlternateSphericalDiffuse<IntegratorPolicy, ProfilePolicy> >(), origin_(o), maxLGreen_(l)
    {
        initProfilePolicy(e1, e2, w, c);
        initSphericalDiffuse();
    }
    virtual ~AlternateSphericalDiffuse() {}
    /*! Returns value of the kernel of the \f$\mathcal{D}\f$ integral operator for the pair of points p1, p2:
     *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this methods with \f$\mathbf{p}_1\f$
     *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} = \mathbf{n}_{\mathbf{p}_1}\f$
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double kernelD(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
    {
        double eps_r2 = 0.0;
        std::tie(eps_r2, std::ignore) = this->epsilon(p2);

        return (eps_r2 * this->derivativeProbe(direction, p1, p2));
    }

    /*! Calculates the matrix representation of the S operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd singleLayer(const std::vector<Element> & e) const override
    {
        return this->integrator_.singleLayer(*this, e);
    }
    /*! Calculates the matrix representation of the D operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd doubleLayer(const std::vector<Element> & e) const override
    {
        return this->integrator_.doubleLayer(*this, e);
    }

    friend std::ostream & operator<<(std::ostream & os, AlternateSphericalDiffuse & gf) {
        return gf.printObject(os);
    }
    /*! \brief Returns Coulomb singularity separation coefficient
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double inverseE(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const {
        // Obtain coefficient for the separation of the Coulomb singularity
        double invE = 0.0;
        for (int L = 1; L < maxLGreen_; ++L) {
            invE += this->inverseE_impl(L, source, probe);
        }
        return invE;
    }
    /*! Returns value of the directional derivative of the
     *  Coulomb singularity separation coefficient for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point, thus assuming that the direction is relative to that point.
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    double inverseEDerivative(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const {
        using namespace std::placeholders;
        return threePointStencil(std::bind(&AlternateSphericalDiffuse<IntegratorPolicy, ProfilePolicy>::inverseE, this, _1, _2),
                p2, p1, direction, this->delta_);
    }
    /*! Handle to the dielectric profile evaluation */
    std::tuple<double, double> epsilon(const Eigen::Vector3d & point) const {
        return this->profile_((point + this->origin_).norm());
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    using RadialSolution = interfaces::RadialSolution;
    /*! Evaluates the Green's function given a pair of points
     *  \param[in] sp the source point
     *  \param[in] pp the probe point
     *
     *  \note This takes care of the origin shift
     */
    virtual Numerical operator()(Numerical * sp, Numerical * pp) const
    {
        // Transfer raw arrays to Eigen vectors using the Map type
        Eigen::Map<Eigen::Matrix<double, 3, 1> > source(sp), probe(pp);

        double invE= 0.0;
        for (int L = 1; L <= maxLGreen_; ++L) {
            invE += this->inverseE_impl(L, source, probe);
        }
        double r12 = (source - probe).norm();

        return invE / r12;
    }
    virtual std::ostream & printObject(std::ostream & os)
    {
        Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, 0, ", ", "\n", "(", ")");
        os << "Green's function type: spherical diffuse" << std::endl;
        os << this->profile_ << std::endl;
        os << "Sphere center        = " << this->origin_.transpose().format(CleanFmt) << std::endl;
        os << "Angular momentum (Green's function)    = " << this->maxLGreen_;
        return os;
    }
    /*! Initializes a one-layer profile
     *  \param[in] e1 left-side dielectric constant
     *  \param[in] e2 right-side dielectric constant
     *  \param[in] w width of the interface layer
     *  \param[in] c center of the diffuse layer
     */
    void initProfilePolicy(double e1, double e2, double w, double c) {
        this->profile_ = ProfilePolicy(e1, e2, w, c);
    }
    /*! This calculates all the components needed to evaluate the Green's function */
    void initSphericalDiffuse() {
        using namespace std::placeholders;
        using namespace interfaces;

        LOG("AlternateSphericalDiffuse::initSphericalDiffuse");
        // Parameters for the numerical solution of the radial differential equation
        double eps_abs_     = 1.0e-10; /*! Absolute tolerance level */
        double eps_rel_     = 1.0e-06; /*! Relative tolerance level */
        double factor_x_    = 0.0;     /*! Weight of the state      */
        double factor_dxdt_ = 0.0;     /*! Weight of the state derivative */
        double r_0_         = 0.5;     /*! Lower bound of the integration interval */
        double r_infinity_  = this->profile_.center() + 200.0; /*! Upper bound of the integration interval */
        double observer_step_ = 1.0e-03; /*! Time step between observer calls */
        IntegratorParameters params_(eps_abs_, eps_rel_, factor_x_, factor_dxdt_, r_0_, r_infinity_, observer_step_);
        ProfileEvaluator eval_ = std::bind(&ProfilePolicy::operator(), this->profile_, _1);

        LOG("Computing radial solutions for Green's function");
        timerON("AlternateSphericalDiffuse: Looping over angular momentum");
        for (int L = 0; L <= maxLGreen_; ++L) {
            // First radial solution
            LOG("Computing first radial solution L = " + std::to_string(L));
            timerON("computeZeta L = " + std::to_string(L));
            // Create an empty RadialSolution
            RadialSolution tmp_zeta_;
            computeZeta(L, tmp_zeta_, eval_, params_);
            zeta_.push_back(tmp_zeta_);
            timerOFF("computeZeta L = " + std::to_string(L));
            LOG("DONE: Computing first radial solution L = " + std::to_string(L));

            // Second radial solution
            LOG("Computing second radial solution L = " + std::to_string(L));
            timerON("computeOmega L = " + std::to_string(L));
            // Create an empty RadialSolution
            RadialSolution tmp_omega_;
            computeOmega(L, tmp_omega_, eval_, params_);
            omega_.push_back(tmp_omega_);
            timerOFF("computeOmega L = " + std::to_string(L));
            LOG("DONE: Computing second radial solution L = " + std::to_string(L));
        }
        timerOFF("AlternateSphericalDiffuse: Looping over angular momentum");
        LOG("DONE: Computing radial solutions for Green's function");
    }

    /*! Center of the dielectric sphere */
    Eigen::Vector3d origin_;

    /**@{ Parameters and functions for the calculation of the Green's function, including Coulomb singularity */
    /*! Maximum angular momentum in the final summation over Legendre polynomials to obtain G */
    int maxLGreen_;
    /*! \brief First independent radial solution, used to build Green's function.
     *  \note The vector has dimension maxLGreen_ and has r^l behavior
     */
    std::vector<RadialSolution> zeta_;
    /*! \brief Second independent radial solution, used to build Green's function.
     *  \note The vector has dimension maxLGreen_  and has r^(-l-1) behavior
     */
    std::vector<RadialSolution> omega_;
    /*! \brief Returns L-th component of the coefficient of the Coulomb singularity
     *  \param[in] L  angular momentum
     *  \param[in] sp source point
     *  \param[in] pp probe point
     *  \note This function shifts the given source and probe points by the location of the
     *  dielectric sphere.
     *
     *  No explicit separation of the Coulomb singularity from the image potential is made.
     *  We rewrite G(r, r') = 1 / [E(r, r') * |r - r'|]
     *  This function returns 1 / E(r, r')
     */
    double inverseE_impl(int L, const Eigen::Vector3d & sp, const Eigen::Vector3d & pp) const {
        using namespace interfaces;
        double r_0_         = 0.5;     /*! Lower bound of the integration interval */
        double r_infinity_  = this->profile_.center() + 200.0; /*! Upper bound of the integration interval */
        Eigen::Vector3d sp_shift = sp + this->origin_;
        Eigen::Vector3d pp_shift = pp + this->origin_;
        double r1 = sp_shift.norm();
        double r2 = pp_shift.norm();
        double cos_gamma = sp_shift.dot(pp_shift) / (r1 * r2);
        double r12 = (sp_shift - pp_shift).norm();
        // Evaluate Legendre polynomial of order L
        // First of all clean-up cos_gamma, Legendre polynomials
        // are only defined for -1 <= x <= 1
        if (numericalZero(cos_gamma - 1)) cos_gamma = 1.0;
        if (numericalZero(cos_gamma + 1)) cos_gamma = -1.0;
        double pl_x = boost::math::legendre_p(L, cos_gamma);

        /* Value of zeta_[L] at point with index 1 */
        double zeta1 = zeta(zeta_[L], L, r1, r_0_);
        /* Value of zeta_[L] at point with index 2 */
        double zeta2 = zeta(zeta_[L], L, r2, r_0_);
        /* Value of omega_[L] at point with index 1 */
        double omega1 = omega(omega_[L], L, r1, r_infinity_);
        /* Value of omega_[L] at point with index 2 */
        double omega2 = omega(omega_[L], L, r2, r_infinity_);

        /* Components for the evaluation of the Wronskian */
        /* Value of derivative of zeta_[L] at point with index 2 */
        double d_zeta2  = derivative_zeta(zeta_[L], L, r2, r_0_);
        /* Value of derivative of omega_[L] at point with index 2 */
        double d_omega2  = derivative_omega(omega_[L], L, r2, r_infinity_);

        double eps_r2 = 0.0;
        std::tie(eps_r2, std::ignore) = this->profile_(pp_shift.norm());

        double denominator = (d_zeta2 - d_omega2) * std::pow(r2, 2) * eps_r2;

        double invE = 0.0;
        if (r1 < r2) {
            invE = std::exp(zeta1 - zeta2) * (2*L +1) * r12 * pl_x / denominator;
        } else {
            invE = std::exp(omega1 - omega2) * (2*L +1) * r12 * pl_x / denominator;
        }

        return invE;
    }
    /**@}*/
};

#endif // ALTERNATESPHERICALDIFFUSE_HPP
