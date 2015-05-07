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

#ifndef SPHERICALDIFFUSE_HPP
#define SPHERICALDIFFUSE_HPP

#include <array>
#include <cmath>
#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>
#include "taylor.hpp"

#include "MathUtils.hpp"

template <typename ProfilePolicy>
class SphericalDiffuse;

/*! \file SphericalDiffuse.hpp
 *  \typedef TanhSphericalDiffuse
 *  \brief Green's function for a diffuse interface with spherical symmetry and tanh profile
 *  \author Hui Cao, Ville Weijo, Luca Frediani and Roberto Di Remigio
 *  \date 2010-2015
 */
typedef SphericalDiffuse<TanhDiffuse> TanhSphericalDiffuse;

#include "DerivativeTypes.hpp"
#include "DiagonalIntegratorFactory.hpp"
#include "DiagonalIntegrator.hpp"
#include "InterfacesDef.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file SphericalDiffuse.hpp
 *  \class SphericalDiffuse
 *  \brief Green's function for a diffuse interface with spherical symmetry
 *  \author Hui Cao, Ville Weijo, Luca Frediani and Roberto Di Remigio
 *  \date 2010-2015
 *  \tparam ProfilePolicy functional form of the diffuse layer
 *
 *  This class is general, in the sense that no specific dielectric
 *  profile has been set in its definition.
 *  In principle any profile that can be described by:
 *  1. a left-side dielectric constant;
 *  2. a right-side dielectric constant;
 *  3. an interface layer width;
 *  4. an interface layer center
 *  can be used to define a new diffuse interface with spherical symmetry.
 */

template <typename ProfilePolicy>
class SphericalDiffuse : public GreensFunction<Numerical, ProfilePolicy>
{
public:
    /*!
     * \param[in] e1 left-side dielectric constant
     * \param[in] e2 right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] c center of the diffuse layer
     */
    SphericalDiffuse(double e1, double e2, double w, double c) : GreensFunction<Numerical, ProfilePolicy>(false)
    {
        initProfilePolicy(e1, e2, w, c);
        initSphericalDiffuse();
    }
    /*!
     * \param[in] e1 left-side dielectric constant
     * \param[in] e2 right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] c center of the diffuse layer
     * \param[in] diag strategy to calculate the diagonal elements of the boundary integral operator
     */
    SphericalDiffuse(double e1, double e2, double w, double c, DiagonalIntegrator * diag)
        : GreensFunction<Numerical, ProfilePolicy>(false, diag)
    {
        initProfilePolicy(e1, e2, w, c);
        initSphericalDiffuse();
    }
    virtual ~SphericalDiffuse() {}
    /*! Returns value of the directional derivative of the
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point, thus assuming that the direction is relative to that point.
     *
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
    {
        /*
        double eps_r2 = 0.0, epsPrime_r2 = 0.0;
        this->profile_(eps_r2, epsPrime_r2, p2.norm());

        return (eps_r2 * this->derivativeProbe(direction, p1, p2));
        */
        using namespace std::placeholders;
        double CoulombDeriv = threePointStencil(std::bind(&SphericalDiffuse::Coulomb, this, _1, _2),
                p2, p1, direction, this->delta_);
        double imageDeriv = threePointStencil(std::bind(&SphericalDiffuse::imagePotential, this, _1, _2),
                p2, p1, direction, this->delta_);

        double eps_r2 = 0.0, epsPrime_r2 = 0.0;
        this->profile_(eps_r2, epsPrime_r2, p2.norm());

        return (eps_r2 * CoulombDeriv + imageDeriv);
    }
    /*!
     *  Calculates the diagonal elements of the S operator: \f$ S_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalS(const Element & e) const
    {
            return this->diagonal_->computeS(this, e);
    }
    /*!
     *  Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalD(const Element & e) const
    {
            return this->diagonal_->computeD(this, e);
    }

    virtual double epsilon() const { return 0.0; }

    friend std::ostream & operator<<(std::ostream & os, SphericalDiffuse & gf) {
        return gf.printObject(os);
    }
    /*! \brief Returns Coulomb singularity separation coefficient
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double coefficientCoulomb(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const;
    /*! \brief Returns singular part of the Green's function
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double Coulomb(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const;
    /*! \brief Returns non-singular part of the Green's function (image potential)
     *  \param[in] source location of the source charge
     *  \param[in] probe location of the probe charge
     */
    double imagePotential(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const;
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
    double coefficientCoulombDerivative(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
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
    double CoulombDerivative(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
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
    double imagePotentialDerivative(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    /*! Handle to the dielectric profile evaluation
     */
    void epsilon(double & v, double & d, double point) const { this->profile_(v, d, point); }
private:
    /*! Evaluates the Green's function given a pair of points
     *
     *  \param[in] sp the source point
     *  \param[in] pp the probe point
     */
    virtual Numerical operator()(Numerical * sp, Numerical * pp) const
    {
        // Transfer raw arrays to Eigen vectors using the Map type
        Eigen::Map<Eigen::Matrix<double, 3, 1> > p1(sp), p2(pp);
        double r1  = p1.norm();
        double r2  = p2.norm();
        double r12 = (p1 - p2).norm();
        double cos_gamma = p1.dot(p2) / (r1 * r2);

        // Obtain coefficient for the separation of the Coulomb singularity
        double Cr12 = this->coefficient(r1, r2);

        double gr12 = 0.0;
        for (int L = 0; L <= maxLGreen_; ++L) {
            gr12 += this->functionSummation(L, r1, r2, cos_gamma, Cr12);
        }

        return (1.0 / (Cr12 * r12) + gr12);
    }
    virtual std::ostream & printObject(std::ostream & os)
    {
        os << "Green's function type: spherical diffuse" << std::endl;
        os << this->profile_ << std::endl;
        os << "Angular momentum (Green's function)    = " << this->maxLGreen_ << std::endl;
        os << "Angular momentum (Coulomb coefficient) = " << this->maxLC_;
        return os;
    }
    /*!
     * \param[in] e1 left-side dielectric constant
     * \param[in] e2 right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] c center of the diffuse layer
     */
    void initProfilePolicy(double e1, double e2, double w, double c)
    { this->profile_ = ProfilePolicy(e1, e2, w, c); }
    /*! This calculates all the components needed to evaluate the Green's function */
    void initSphericalDiffuse();

    /**@{ Parameters and functions for the calculation of the Green's function, including Coulomb singularity */
    /*! Maximum angular momentum in the final summation over Legendre polynomials to obtain G */
    int maxLGreen_ = 30;
    /*! \brief First independent radial solution, used to build Green's function.
     *  \note The vector has dimension maxLGreen_ and has r^l behavior
     */
    std::vector<RadialFunction> zeta_;
    /*! \brief Second independent radial solution, used to build Green's function.
     *  \note The vector has dimension maxLGreen_  and has r^(-l-1) behavior
     */
    std::vector<RadialFunction> omega_;
    /*! \brief Returns L-th component of the radial part of the Green's function
     *  \param[in] L  angular momentum
     *  \param[in] r1 first point
     *  \param[in] r2 second point
     *  \param[in] cos_gamma cosine of the angle in between first and second point
     *  \param[in] Cr12 Coulomb singularity separation coefficient
     */
    double functionSummation(int L, double r1, double r2, double cos_gamma, double Cr12) const;
    /**@}*/

    /**@{ Parameters and functions for the calculation of the Coulomb singularity separation coefficient */
    /*! Maximum angular momentum to obtain C(r, r'), needed to separate the Coulomb singularity */
    int maxLC_     = 2 * maxLGreen_;
    /*! \brief First independent radial solution, used to build coefficient.
     *  \note This is needed to separate the Coulomb singularity and has r^l behavior
     */
    RadialFunction zetaC_;
    /*! \brief Second independent radial solution, used to build coefficient.
     *  \note This is needed to separate the Coulomb singularity and has r^(-l-1) behavior
     */
    RadialFunction omegaC_;
    /*! \brief Returns coefficient for the separation of the Coulomb singularity
     *  \param[in] r1 first point
     *  \param[in] r2 second point
     */
    double coefficient(double r1, double r2) const;
    /**@}*/

    // This is needed for the semi-analytic implementation
    //Eigen::Vector3d coefficientGradient(const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    //Eigen::Vector3d functionSummationGradient(int L, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2, const Eigen::Vector3d & Cr12_grad) const;
};

#endif // SPHERICALDIFFUSE_HPP
