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

#include <cmath>
#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>
#include "taylor.hpp"

template <typename DerivativeTraits>
class IonicLiquid;

#include "DerivativeTypes.hpp"
#include "DiagonalIntegratorFactory.hpp"
#include "DiagonalIntegrator.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file IonicLiquid.hpp
 *  \class IonicLiquid
 *  \brief Green's functions for ionic liquid, described by the linearized Poisson-Boltzmann equation.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013-2014
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 */

template <typename DerivativeTraits>
class IonicLiquid : public GreensFunction<DerivativeTraits, Yukawa>
{
public:
    IonicLiquid(double epsilon, double kappa) : GreensFunction<DerivativeTraits, Yukawa>(false) { initProfilePolicy(epsilon, kappa); }
    IonicLiquid(double epsilon, double kappa, DiagonalIntegrator * diag) : GreensFunction<DerivativeTraits, Yukawa>(false, diag) { initProfilePolicy(epsilon, kappa); }
    virtual ~IonicLiquid() {}
    /*!
     *  Returns value of the directional derivative of the
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
        return this->profile_.epsilon * (this->derivativeProbe(direction, p1, p2));
    }

    /*!
     *  Calculates the diagonal elements of the S operator: \f$ S_{ii} \f$
     *  \param[in] area   area of the i-th tessera to be calculated
     */
    virtual double diagonalS(double area) const
    {
            this->diagonal_->computeS(this, area);
            return 1.0;
    }
    /*!
     *  Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *  \param[in] area   area of the i-th tessera to be calculated
     *  \param[in] radius radius of the sphere the tessera belongs to
     */
    virtual double diagonalD(double area, double radius) const
    {
            this->diagonal_->computeD(this, area, radius);
            return 1.0;
    }

    virtual double epsilon() const { return this->profile_.epsilon; }

    friend std::ostream & operator<<(std::ostream & os, IonicLiquid & gf) {
        return gf.printObject(os);
    }
private:
    /*!
     *  Evaluates the Green's function given a pair of points
     *
     *  \param[in] sp the source point
     *  \param[in] pp the probe point
     */
    virtual DerivativeTraits operator()(DerivativeTraits * sp, DerivativeTraits * pp) const
    {
	double eps = this->profile_.epsilon;
	double k = this->profile_.kappa; 
        DerivativeTraits distance = sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) +
                          (sp[1] - pp[1]) * (sp[1] - pp[1]) +
                          (sp[2] - pp[2]) * (sp[2] - pp[2]));
        return (exp(-k * distance) / (eps * distance));
    }
    virtual std::ostream & printObject(std::ostream & os)
    {
        os << "Green's function type: ionic liquid" << std::endl;
        os << "Permittivity         = " << this->profile_.epsilon << std::endl;
        os << "Inverse Debye length = " << this->profile_.kappa;
        return os;
    }
    void initProfilePolicy(double eps, double k) { this->profile_ = Yukawa(eps, k); }
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

template <>
inline double GreensFunction<double, Yukawa>::function(const Eigen::Vector3d & source,
                                        const Eigen::Vector3d & probe) const
{
    double sp[3], pp[3], res;
    sp[0] = source(0); sp[1] = source(1); sp[2] = source(2);
    pp[0] = probe(0);  pp[1] = probe(1);  pp[2] = probe(2);
    res = this->operator()(sp, pp);
    return res;
}

template <>
inline double GreensFunction<double, Yukawa>::derivativeSource(const Eigen::Vector3d & normal_p1,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d deltaPlus  = p1 + normal_p1 * delta_ / normal_p1.norm();
    Eigen::Vector3d deltaMinus = p1 - normal_p1 * delta_ / normal_p1.norm();
    double funcPlus  = function(deltaPlus,  p2);
    double funcMinus = function(deltaMinus, p2);
    return (funcPlus - funcMinus)/(2.0*delta_);
}

template <>
inline double GreensFunction<double, Yukawa>::derivativeProbe(const Eigen::Vector3d & normal_p2,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d deltaPlus  = p2 + normal_p2 * delta_ / normal_p2.norm();
    Eigen::Vector3d deltaMinus = p2 - normal_p2 * delta_ / normal_p2.norm();
    double funcPlus  = function(p1, deltaPlus);
    double funcMinus = function(p1, deltaMinus);
    return (funcPlus - funcMinus)/(2.0*delta_);
}

#endif // IONICLIQUID_HPP
