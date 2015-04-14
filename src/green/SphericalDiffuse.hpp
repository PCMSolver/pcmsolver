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

#include <cmath>
#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>
#include "taylor.hpp"

template <typename ProfilePolicy>
class SphericalDiffuse;

#include "DerivativeTypes.hpp"
#include "DiagonalIntegratorFactory.hpp"
#include "DiagonalIntegrator.hpp"
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
     * \param[in] eL left-side dielectric constant
     * \param[in] eR right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] c center of the diffuse layer
     */
    SphericalDiffuse(double eL, double eR, double w, double c) : GreensFunction<Numerical, ProfilePolicy>(false)
    {
        initProfilePolicy(eL, eR, w, c);
    }
    /*!
     * \param[in] eL left-side dielectric constant
     * \param[in] eR right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] c center of the diffuse layer
     */
    SphericalDiffuse(double eL, double eR, double w, double c,
            DiagonalIntegrator * diag) : GreensFunction<Numerical, ProfilePolicy>(false, diag)
    {
        initProfilePolicy(eL, eR, w, c);
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
        return this->derivativeProbe(direction, p1, p2);
    }

    /*! Calculates the diagonal elements of the S operator: \f$ S_{ii} \f$
     *
     *  \param[in] area   area of the i-th tessera
     */
    virtual double diagonalS(double area) const
    {
            return this->diagonal_->computeS(this, area);
    }
    /*! Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *
     *  \param[in] area   area of the i-th tessera
     *  \param[in] radius radius of the sphere the tessera belongs to
     */
    virtual double diagonalD(double area, double radius) const
    {
            return this->diagonal_->computeD(this, area, radius);
    }

    virtual double epsilon() const { return this->profile_.epsilon; }

    friend std::ostream & operator<<(std::ostream & os, SphericalDiffuse & gf) {
        return gf.printObject(os);
    }
private:
    /*! Evaluates the Green's function given a pair of points
     *
     *  \param[in] sp the source point
     *  \param[in] pp the probe point
     */
    virtual Numerical operator()(Numerical * sp, Numerical * pp) const
    {
        Numerical res = 0.0;
        return res;
    }
    virtual std::ostream & printObject(std::ostream & os)
    {
        os << "Green's function type: spherical diffuse";
        return os;
    }
    /*!
     * \param[in] eL left-side dielectric constant
     * \param[in] eR right-side dielectric constant
     * \param[in] w width of the interface layer
     * \param[in] c center of the diffuse layer
     */
    void initProfilePolicy(double eL, double eR, double w, double c)
    { this->profile_ = ProfilePolicy(eL, eR, w, c); }
};

#endif // SPHERICALDIFFUSE_HPP
