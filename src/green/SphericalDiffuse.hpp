/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#ifndef SPHERICALDIFFUSE_HPP
#define SPHERICALDIFFUSE_HPP

#include <cmath>
#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

// Has to be included here
#include "InterfacesImpl.hpp"
// Boost.Math includes
#include <boost/math/special_functions/legendre.hpp>

#include "DerivativeUtils.hpp"
#include "bi_operators/IntegratorForward.hpp"
#include "dielectric_profile/ProfileForward.hpp"
#include "GreensFunction.hpp"
#include "utils/MathUtils.hpp"

/*! \file SphericalDiffuse.hpp
 *  \class SphericalDiffuse
 *  \brief Green's function for a diffuse interface with spherical symmetry
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

template <typename IntegratorPolicy = CollocationIntegrator,
          typename ProfilePolicy = OneLayerTanh>
class SphericalDiffuse __final : public GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy,
                                               SphericalDiffuse<IntegratorPolicy, ProfilePolicy> >
{
public:
    /*! Constructor for a one-layer interface
     * \param[in] e1 left-side dielectric constant
     * \param[in] e2 right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] c center of the diffuse layer
     * \param[in] o center of the sphere
     */
    SphericalDiffuse(double e1, double e2, double w, double c, const Eigen::Vector3d & o, int l)
        : GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy, SphericalDiffuse<IntegratorPolicy, ProfilePolicy> >(),
          origin_(o), maxLGreen_(l), maxLC_(2*l)
    {
        initProfilePolicy(e1, e2, w, c);
        initSphericalDiffuse();
    }
    /*! Constructor for a one-layer interface
     * \param[in] e1 left-side dielectric constant
     * \param[in] e2 right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] c center of the diffuse layer
     * \param[in] o center of the sphere
     */
    SphericalDiffuse(double e1, double e2, double w, double c, const Eigen::Vector3d & o, int l, double f)
        : GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy, SphericalDiffuse<IntegratorPolicy, ProfilePolicy> >(f),
          origin_(o), maxLGreen_(l), maxLC_(2*l)
    {
        initProfilePolicy(e1, e2, w, c);
        initSphericalDiffuse();
    }
    virtual ~SphericalDiffuse() {}

    /*! Calculates the matrix representation of the S operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd singleLayer(const std::vector<Element> & e) const __override
    {
        return this->integrator_.singleLayer(*this, e);
    }
    /*! Calculates the matrix representation of the D operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd doubleLayer(const std::vector<Element> & e) const __override
    {
        return this->integrator_.doubleLayer(*this, e);
    }

    friend std::ostream & operator<<(std::ostream & os, SphericalDiffuse & gf) {
        return gf.printObject(os);
    }
    /*! \brief Returns Coulomb singularity separation coefficient
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double coefficientCoulomb(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const {
        // Obtain coefficient for the separation of the Coulomb singularity
        return this->coefficient_impl(source, probe);
    }
    /*! \brief Returns singular part of the Green's function
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double Coulomb(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const {
        double r12 = (source - probe).norm();

        // Obtain coefficient for the separation of the Coulomb singularity
        return (1.0 / (this->coefficient_impl(source, probe) * r12));
    }
    /*! \brief Returns non-singular part of the Green's function (image potential)
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double imagePotential(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const {
        // Obtain coefficient for the separation of the Coulomb singularity
        double Cr12 = this->coefficient_impl(source, probe);

        double gr12 = 0.0;
        for (int L = 1; L <= maxLGreen_; ++L) {
            gr12 += this->imagePotentialComponent_impl(L, source, probe, Cr12);
        }

        return gr12;
    }
    /*! Returns value of the directional derivative of the
     *  Coulomb singularity separation coefficient for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point, thus assuming that the direction is relative to that point.
     *
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    double coefficientCoulombDerivative(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const {
        return threePointStencil(pcm::bind(&SphericalDiffuse<IntegratorPolicy, ProfilePolicy>::coefficientCoulomb, this, pcm::_1, pcm::_2),
                p2, p1, direction, this->delta_);
    }
    /*! Returns value of the directional derivative of the
     *  singular part of the Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point, thus assuming that the direction is relative to that point.
     *
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    double CoulombDerivative(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const {
        return threePointStencil(pcm::bind(&SphericalDiffuse<IntegratorPolicy, ProfilePolicy>::Coulomb, this, pcm::_1, pcm::_2),
                p2, p1, direction, this->delta_);
    }
    /*! Returns value of the directional derivative of the
     *  non-singular part (image potential) of the Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point, thus assuming that the direction is relative to that point.
     *
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    double imagePotentialDerivative(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const {
        return threePointStencil(pcm::bind(&SphericalDiffuse<IntegratorPolicy, ProfilePolicy>::imagePotential, this, pcm::_1, pcm::_2),
                p2, p1, direction, this->delta_);
    }
    /*! Handle to the dielectric profile evaluation */
    pcm::tuple<double, double> epsilon(const Eigen::Vector3d & point) const {
        return this->profile_((point + this->origin_).norm());
    }
    void toFile(const std::string & prefix = "") {
	std::string tmp;
	prefix.empty() ? tmp = prefix : tmp = prefix + "-";
        writeToFile(zetaC_, tmp + "zetaC.dat");
        writeToFile(omegaC_, tmp + "omegaC.dat");
	for (int L = 1; L <= maxLGreen_; ++L) {
	    writeToFile(zeta_[L], tmp + "zeta_" + pcm::to_string(L) + ".dat");
	    writeToFile(omega_[L], tmp + "omega_" + pcm::to_string(L) + ".dat");
	}
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    /*! Evaluates the Green's function given a pair of points
     *  \param[in] sp the source point
     *  \param[in] pp the probe point
     *
     *  \note This takes care of the origin shift
     */
    virtual Numerical operator()(Numerical * sp, Numerical * pp) const __override
    {
        // Transfer raw arrays to Eigen vectors using the Map type
        Eigen::Map<Eigen::Matrix<double, 3, 1> > source(sp), probe(pp);

        // Obtain coefficient for the separation of the Coulomb singularity
        double Cr12 = this->coefficient_impl(source, probe);

        double gr12 = 0.0;
        for (int L = 1; L <= maxLGreen_; ++L) {
            gr12 += this->imagePotentialComponent_impl(L, source, probe, Cr12);
        }
        double r12 = (source - probe).norm();

        return (1.0 / (Cr12 * r12) + gr12);
    }
    /*! Returns value of the kernel of the \f$\mathcal{D}\f$ integral operator for the pair of points p1, p2:
     *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this methods with \f$\mathbf{p}_1\f$
     *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} = \mathbf{n}_{\mathbf{p}_1}\f$
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double kernelD_impl(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const __override
    {
        double eps_r2 = 0.0;
        // Shift p2 by origin_
        pcm::tie(eps_r2, pcm::ignore) = this->epsilon(p2);

        return (eps_r2 * this->derivativeProbe(direction, p1, p2));
    }
    virtual KernelS exportKernelS_impl() const __override {
      return pcm::bind(&SphericalDiffuse<IntegratorPolicy, ProfilePolicy>::kernelS, *this, pcm::_1, pcm::_2);
    }
    virtual KernelD exportKernelD_impl() const __override {
      return pcm::bind(&SphericalDiffuse<IntegratorPolicy, ProfilePolicy>::kernelD, *this, pcm::_1, pcm::_2, pcm::_3);
    }
    virtual std::ostream & printObject(std::ostream & os) __override
    {
        Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, 0, ", ", "\n", "(", ")");
        os << "Green's function type: spherical diffuse" << std::endl;
        os << this->profile_ << std::endl;
        os << "Sphere center        = " << this->origin_.transpose().format(CleanFmt) << std::endl;
        os << "Angular momentum (Green's function)    = " << this->maxLGreen_ << std::endl;
        os << "Angular momentum (Coulomb coefficient) = " << this->maxLC_;
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
        using namespace interfaces;

        LOG("SphericalDiffuse::initSphericalDiffuse");
        // Parameters for the numerical solution of the radial differential equation
        double eps_abs_     = 1.0e-10; /*! Absolute tolerance level */
        double eps_rel_     = 1.0e-06; /*! Relative tolerance level */
        double factor_x_    = 0.0;     /*! Weight of the state      */
        double factor_dxdt_ = 0.0;     /*! Weight of the state derivative */
        double r_0_         = 0.5;     /*! Lower bound of the integration interval */
        double r_infinity_  = this->profile_.center() + 200.0; /*! Upper bound of the integration interval */
        double observer_step_ = 1.0e-03; /*! Time step between observer calls */
        IntegratorParameters params_(eps_abs_, eps_rel_, factor_x_, factor_dxdt_, r_0_, r_infinity_, observer_step_);
        ProfileEvaluator eval_ = pcm::bind(&ProfilePolicy::operator(), this->profile_, pcm::_1);

        LOG("Computing coefficient for the separation of the Coulomb singularity");
        LOG("Computing first radial solution L = " + pcm::to_string(maxLC_));
        TIMER_ON("computeZeta for coefficient");
        zetaC_ = RadialFunction<StateType, LnTransformedRadial, Zeta>(maxLC_, r_0_, r_infinity_, eval_, params_);
        TIMER_OFF("computeZeta for coefficient");
        LOG("DONE: Computing first radial solution L = " + pcm::to_string(maxLC_));

        LOG("Computing second radial solution L = " + pcm::to_string(maxLC_));
        TIMER_ON("computeOmega for coefficient");
        omegaC_ = RadialFunction<StateType, LnTransformedRadial, Omega>(maxLC_, r_0_, r_infinity_, eval_, params_);
        TIMER_OFF("computeOmega for coefficient");
        LOG("Computing second radial solution L = " + pcm::to_string(maxLC_));
        LOG("DONE: Computing coefficient for the separation of the Coulomb singularity");

        LOG("Computing radial solutions for Green's function");
        TIMER_ON("SphericalDiffuse: Looping over angular momentum");
        zeta_.reserve(maxLGreen_+1);
        omega_.reserve(maxLGreen_+1);
        for (int L = 0; L <= maxLGreen_; ++L) {
            // First radial solution
            LOG("Computing first radial solution L = " + pcm::to_string(L));
            TIMER_ON("computeZeta L = " + pcm::to_string(L));
            // Create an empty RadialSolution
            RadialFunction<StateType, LnTransformedRadial, Zeta> tmp_zeta_(L, r_0_, r_infinity_, eval_, params_);
            zeta_.push_back(tmp_zeta_);
            TIMER_OFF("computeZeta L = " + pcm::to_string(L));
            LOG("DONE: Computing first radial solution L = " + pcm::to_string(L));

            // Second radial solution
            LOG("Computing second radial solution L = " + pcm::to_string(L));
            TIMER_ON("computeOmega L = " + pcm::to_string(L));
            // Create an empty RadialSolution
            RadialFunction<StateType, LnTransformedRadial, Omega> tmp_omega_(L, r_0_, r_infinity_, eval_, params_);
            omega_.push_back(tmp_omega_);
            TIMER_OFF("computeOmega L = " + pcm::to_string(L));
            LOG("DONE: Computing second radial solution L = " + pcm::to_string(L));
        }
        TIMER_OFF("SphericalDiffuse: Looping over angular momentum");
        LOG("DONE: Computing radial solutions for Green's function");
    }

    /*! Center of the dielectric sphere */
    Eigen::Vector3d origin_;

    /**@{ Parameters and functions for the calculation of the Green's function, including Coulomb singularity */
    /*! Maximum angular momentum in the __final summation over Legendre polynomials to obtain G */
    int maxLGreen_;
    /*! \brief First independent radial solution, used to build Green's function.
     *  \note The vector has dimension maxLGreen_ and has r^l behavior
     */
    std::vector<RadialFunction<interfaces::StateType, interfaces::LnTransformedRadial, Zeta> > zeta_;
    /*! \brief Second independent radial solution, used to build Green's function.
     *  \note The vector has dimension maxLGreen_  and has r^(-l-1) behavior
     */
    std::vector<RadialFunction<interfaces::StateType, interfaces::LnTransformedRadial, Omega> > omega_;
    /*! \brief Returns L-th component of the radial part of the Green's function
     *  \param[in] L  angular momentum
     *  \param[in] sp source point
     *  \param[in] pp probe point
     *  \param[in] Cr12 Coulomb singularity separation coefficient
     *  \note This function shifts the given source and probe points by the location of the
     *  dielectric sphere.
     */
    double imagePotentialComponent_impl(int L, const Eigen::Vector3d & sp, const Eigen::Vector3d & pp, double Cr12) const {
        Eigen::Vector3d sp_shift = sp + this->origin_;
        Eigen::Vector3d pp_shift = pp + this->origin_;
        double r1 = sp_shift.norm();
        double r2 = pp_shift.norm();
        double cos_gamma = sp_shift.dot(pp_shift) / (r1 * r2);
        // Evaluate Legendre polynomial of order L
        // First of all clean-up cos_gamma, Legendre polynomials
        // are only defined for -1 <= x <= 1
        if (numericalZero(cos_gamma - 1)) cos_gamma = 1.0;
        if (numericalZero(cos_gamma + 1)) cos_gamma = -1.0;
        double pl_x = boost::math::legendre_p(L, cos_gamma);

        /* Sample zeta_[L] */
        double zeta1 = 0.0, zeta2 = 0.0, d_zeta2 = 0.0;
        /* Value of zeta_[L] at point with index 1 */
        pcm::tie(zeta1, pcm::ignore) = zeta_[L](r1);
        /* Value of zeta_[L] and its first derivative at point with index 2 */
        pcm::tie(zeta2, d_zeta2) = zeta_[L](r2);

        /* Sample omega_[L] */
        double omega1 = 0.0, omega2 = 0.0, d_omega2 = 0.0;
        /* Value of omega_[L] at point with index 1 */
        pcm::tie(omega1, pcm::ignore) = omega_[L](r1);
        /* Value of omega_[L] and its first derivative at point with index 2 */
        pcm::tie(omega2, d_omega2) = omega_[L](r2);

        double eps_r2 = 0.0;
        pcm::tie(eps_r2, pcm::ignore) = this->profile_(pp_shift.norm());

        /* Evaluation of the Wronskian and the denominator */
        double denominator = (d_zeta2 - d_omega2) * std::pow(r2, 2) * eps_r2;

        double gr12 = 0.0;
        if (r1 < r2) {
            gr12 = std::exp(zeta1 - zeta2) * (2*L +1) / denominator;
	    double f_L = r1 / r2;
	    for (int i = 1; i < L; ++i) { f_L *= r1 / r2; }
            gr12 = (gr12 - f_L / (r2 * Cr12) ) * pl_x ;
        } else {
            gr12 = std::exp(omega1 - omega2) * (2*L +1) / denominator;
	    double f_L = r2 / r1;
	    for (int i = 1; i < L; ++i) { f_L *= r2 / r1; }
            gr12 = (gr12 - f_L / (r1 * Cr12) ) * pl_x ;
        }

        return gr12;
    }
    /**@}*/

    /**@{ Parameters and functions for the calculation of the Coulomb singularity separation coefficient */
    /*! Maximum angular momentum to obtain C(r, r'), needed to separate the Coulomb singularity */
    int maxLC_; // = 2 * maxLGreen_;
    /*! \brief First independent radial solution, used to build coefficient.
     *  \note This is needed to separate the Coulomb singularity and has r^l behavior
     */
    RadialFunction<interfaces::StateType, interfaces::LnTransformedRadial, Zeta> zetaC_;
    /*! \brief Second independent radial solution, used to build coefficient.
     *  \note This is needed to separate the Coulomb singularity and has r^(-l-1) behavior
     */
    RadialFunction<interfaces::StateType, interfaces::LnTransformedRadial, Omega> omegaC_;
    /*! \brief Returns coefficient for the separation of the Coulomb singularity
     *  \param[in] sp first point
     *  \param[in] pp second point
     *  \note This function shifts the given source and probe points by the location of the
     *  dielectric sphere.
     */
    double coefficient_impl(const Eigen::Vector3d & sp, const Eigen::Vector3d & pp) const {
        double r1 = (sp + this->origin_).norm();
        double r2 = (pp + this->origin_).norm();

        /* Sample zetaC_ */
        double zeta1 = 0.0, zeta2 = 0.0, d_zeta2 = 0.0;
        /* Value of zetaC_ at point with index 1 */
        pcm::tie(zeta1, pcm::ignore) = zetaC_(r1);
        /* Value of zetaC_ and its first derivative at point with index 2 */
        pcm::tie(zeta2, d_zeta2) = zetaC_(r2);

        /* Sample omegaC_ */
        double omega1 = 0.0, omega2 = 0.0, d_omega2 = 0.0;
        /* Value of omegaC_ at point with index 1 */
        pcm::tie(omega1, pcm::ignore) = omegaC_(r1);
        /* Value of omegaC_ and its first derivative at point with index 2 */
        pcm::tie(omega2, d_omega2) = omegaC_(r2);

        double tmp = 0.0, coeff = 0.0;
        double eps_r2 = 0.0;
        pcm::tie(eps_r2, pcm::ignore) = this->profile_(r2);

        /* Evaluation of the Wronskian and the denominator */
        double denominator = (d_zeta2 - d_omega2) * std::pow(r2, 2) * eps_r2;

        if (r1 < r2) {
	    double f_L = r1 / r2;
	    for (int i = 1; i < maxLC_; ++i) { f_L *= r1 / r2; }
            tmp = std::exp(zeta1 - zeta2) * (2*maxLC_ +1) / denominator;
            coeff = f_L / (tmp * r2);
        } else {
	    double f_L = r2 / r1;
	    for (int i = 1; i < maxLC_; ++i) { f_L *= r2 / r1; }
            tmp = std::exp(omega1 - omega2) * (2*maxLC_ +1) / denominator;
            coeff = f_L / (tmp * r1);
        }

        return coeff;
    }
    /**@}*/
};

#endif // SPHERICALDIFFUSE_HPP
