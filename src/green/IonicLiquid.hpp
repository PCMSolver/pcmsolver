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

#ifndef IONICLIQUID_HPP
#define IONICLIQUID_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

class DiagonalIntegrator;

#include "DerivativeTypes.hpp"
#include "DiagonalIntegratorFactory.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file IonicLiquid.hpp
 *  \class IonicLiquid
 *  \brief Green's functions for ionic liquid, described by the linearized Poisson-Boltzmann equation.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013-2014
 *  \tparam T evaluation strategy for the function and its derivatives
 */

template <typename T>
class IonicLiquid : public GreensFunction<T>
{
public:
    IonicLiquid(double epsilon, double kappa) : GreensFunction<T>(false),
        epsilon_(epsilon), kappa_(kappa) {}
    IonicLiquid(double epsilon, double kappa, DiagonalIntegrator * diag) : GreensFunction<T>(false, diag),
        epsilon_(epsilon), kappa_(kappa) {}
    virtual ~IonicLiquid() {}
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

    virtual double epsilon() const { return epsilon_; }

    friend std::ostream & operator<<(std::ostream & os, IonicLiquid & gf) {
        return gf.printObject(os);
    }
private:
    /*!
     *  Evaluates the Green's function given a pair of points
     *
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual T operator()(T source[3], T probe[3]) const;
    double epsilon_;
    double kappa_;
    virtual std::ostream & printObject(std::ostream & os);
};

namespace
{
    struct buildIonicLiquid {
        template <typename DerivativeType>
        IGreensFunction * operator()(const greenData & _data) {
            DiagonalIntegrator * integrator =
		    DiagonalIntegratorFactory::TheDiagonalIntegratorFactory().createDiagonalIntegrator(_data.integratorType);
            return new IonicLiquid<DerivativeType>(_data.epsilon, _data.kappa, integrator);
        }
    };

    IGreensFunction * createIonicLiquid(const greenData & _data)
    {
        buildIonicLiquid build;
        return for_id<derivative_types>(build, _data, _data.how);
    }
    const std::string IONICLIQUID("IONICLIQUID");
    const bool registeredIonicLiquid =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            IONICLIQUID, createIonicLiquid);
}

#endif // IONICLIQUID_HPP
