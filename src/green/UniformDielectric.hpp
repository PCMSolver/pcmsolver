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

#ifndef UNIFORMDIELECTRIC_HPP
#define UNIFORMDIELECTRIC_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

class DiagonalIntegrator;

#include "DerivativeTypes.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file UniformDielectric.hpp
 *  \class UniformDielectric
 *  \brief Green's function for uniform dielectric.
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2014
 *  \tparam T evaluation strategy for the function and its derivatives
 */

template <typename T>
class UniformDielectric : public GreensFunction<T>
{
public:
    UniformDielectric(double eps) : GreensFunction<T>(true), epsilon_(eps) {}
    explicit UniformDielectric(double eps, DiagonalIntegrator * diag) : GreensFunction<T>(true, diag), epsilon_(eps) {}
    virtual ~UniformDielectric() {}
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
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;

    /*!
     *  Calculates the diagonal elements of the S operator: \f$ S_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalS(const Element & e) const;
    /*!
     *  Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalD(const Element & e) const;

    virtual void epsilon(double eps) { epsilon_ = eps; }
    virtual double epsilon() const { return epsilon_; }

    friend std::ostream & operator<<(std::ostream & os, UniformDielectric & gf) {
        return gf.printObject(os);
    }
private:
    /*!
     *  Evaluates the Green's function given a pair of points
     *
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual T operator()(T * source, T * probe) const;
    virtual std::ostream & printObject(std::ostream & os);
    double epsilon_;
};

namespace
{
    struct buildUniformDielectric {
        template <typename DerivativeType>
        IGreensFunction * operator()(const greenData & _data) {
            return new UniformDielectric<DerivativeType>(_data.epsilon, _data.integrator);
        }
    };

    IGreensFunction * createUniformDielectric(const greenData & _data)
    {
        buildUniformDielectric build;
        return for_id<derivative_types>(build, _data, _data.how);
    }
    const std::string UNIFORMDIELECTRIC("UNIFORMDIELECTRIC");
    const bool registeredUniformDielectric =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            UNIFORMDIELECTRIC, createUniformDielectric);
}

#endif // UNIFORMDIELECTRIC_HPP
