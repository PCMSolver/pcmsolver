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

#ifndef GREENSFUNCTION_HPP
#define GREENSFUNCTION_HPP

#include <cmath>
#include <functional>
#include <iosfwd>
#include <stdexcept>

#include "Config.hpp"

#include <Eigen/Dense>

class DiagonalIntegrator;
class Element;

#include "DerivativeTypes.hpp"
#include "IGreensFunction.hpp"
#include "MathUtils.hpp"
#include "ProfileTypes.hpp"

/*! \file GreensFunction.hpp
 *  \class GreensFunction
 *  \brief Templated interface for Green's functions
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2014
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam ProfilePolicy    dielectric profile type
 */

template <typename DerivativeTraits = AD_directional,
	  typename ProfilePolicy = Uniform>
class GreensFunction: public IGreensFunction
{
public:
    GreensFunction(bool uniform) : IGreensFunction(uniform), delta_(1.0e-4) {}
    GreensFunction(bool uniform, DiagonalIntegrator * diag) : IGreensFunction(uniform, diag), delta_(1.0e-4) {}
    virtual ~GreensFunction() {}
    /*!
     *  Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e. the value of the
     *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    virtual double function(const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
    {
        DerivativeTraits sp[3], pp[3], res;
        sp[0] = p1(0); sp[1] = p1(1); sp[2] = p1(2);
        pp[0] = p2(0); pp[1] = p2(1); pp[2] = p2(2);
        res = this->operator()(sp, pp);
        return res[0];
    }
    /*!
     *  Returns value of the kernel for the calculation of the \f$\mathcal{D}\f$ integral operator
     *  for the pair of points p1, p2:
     *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this methods with \f$\mathbf{p}_1\f$
     *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} = \mathbf{n}_{\mathbf{p}_1}\f$
     *
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const = 0;
    /*! Returns value of the directional derivative of the
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_1}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the source point.
     *
     *  \param[in] normal_p1 the normal vector to p1
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivativeSource(const Eigen::Vector3d & normal_p1,
                            const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
    {
        DerivativeTraits t1[3], t2[3], der;
        t1[0] = p1(0); t1[1] = p1(1); t1[2] = p1(2);
        t1[0][1] = normal_p1(0); t1[1][1] = normal_p1(1); t1[2][1] = normal_p1(2);
        t2[0] = p2(0); t2[1] = p2(1); t2[2] = p2(2);
        der = this->operator()(t1, t2);
        return der[1];
    }
    /*! Returns value of the directional derivative of the
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point.
     *
     *  \param[in] normal_p2 the normal vector to p2
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivativeProbe(const Eigen::Vector3d & normal_p2,
                                   const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
    {
        DerivativeTraits t1[3], t2[3], der;
        t1[0] = p1(0); t1[1] = p1(1); t1[2] = p1(2);
        t2[0] = p2(0); t2[1] = p2(1); t2[2] = p2(2);
        t2[0][1] = normal_p2(0); t2[1][1] = normal_p2(1); t2[2][1] = normal_p2(2);
        der = this->operator()(t1, t2);
        return der[1];
    }

    /*!
     *  Calculates the diagonal elements of the S operator: \f$ S_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalS(const Element & e) const = 0;
    /*!
     *  Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalD(const Element & e) const = 0;

    virtual void delta(double value)
    {
        if (value <= 1.0e-10) {
            throw std::invalid_argument("Delta value must be larger than 1.0e-10");
        }
        delta_ = value;
    }
    virtual double delta() { return delta_; }
    virtual double epsilon() const = 0;

    friend std::ostream & operator<<(std::ostream & os, GreensFunction & gf) {
        return gf.printObject(os);
    }
protected:
    /*! Evaluates the Green's function given a pair of points
     *
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual DerivativeTraits operator()(DerivativeTraits * source, DerivativeTraits * probe) const = 0;
    virtual std::ostream & printObject(std::ostream & os)
    {
        os << "Green's Function" << std::endl;
        os << "Delta = " << delta_ << std::endl;
        os << "Uniform = " << uniform_;
        return os;
    }
    double delta_;
    ProfilePolicy profile_;
};

template <>
inline double GreensFunction<Numerical, Uniform>::function(const Eigen::Vector3d & source,
                                        const Eigen::Vector3d & probe) const
{
    Numerical sp[3], pp[3], res;
    sp[0] = source(0); sp[1] = source(1); sp[2] = source(2);
    pp[0] = probe(0);  pp[1] = probe(1);  pp[2] = probe(2);
    res = this->operator()(sp, pp);
    return res;
}

template <>
inline double GreensFunction<Numerical, Uniform>::derivativeSource(const Eigen::Vector3d & normal_p1,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    using namespace std::placeholders;
    return threePointStencil(std::bind(&GreensFunction<Numerical, Uniform>::function, this, _1, _2),
                            p1, p2, normal_p1, this->delta_);
}

template <>
inline double GreensFunction<Numerical, Uniform>::derivativeProbe(const Eigen::Vector3d & normal_p2,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    using namespace std::placeholders;
    return threePointStencil(std::bind(&GreensFunction<Numerical, Uniform>::function, this, _1, _2),
                            p2, p1, normal_p2, this->delta_);
}

#endif // GREENSFUNCTION_HPP
