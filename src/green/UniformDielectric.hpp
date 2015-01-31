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

#include <cmath>
#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>
#include "taylor.hpp"

template <typename DerivativeTraits>
class UniformDielectric;

#include "DerivativeTypes.hpp"
#include "DiagonalIntegratorFactory.hpp"
#include "DiagonalIntegrator.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file UniformDielectric.hpp
 *  \class UniformDielectric
 *  \brief Green's function for uniform dielectric.
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2014
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 */

template <typename DerivativeTraits>
class UniformDielectric : public GreensFunction<DerivativeTraits>
{
public:
    UniformDielectric(double eps) : GreensFunction<DerivativeTraits>(true), epsilon_(eps) {}
    explicit UniformDielectric(double eps, DiagonalIntegrator * diag) : GreensFunction<DerivativeTraits>(true, diag), epsilon_(eps) {}
    virtual ~UniformDielectric() {}
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
        return epsilon_ * (this->derivativeProbe(direction, p1, p2));
    }

    /*! Calculates the diagonal elements of the S operator: \f$ S_{ii} \f$
     *  \param[in] area   area of the i-th tessera to be calculated
     */
    virtual double diagonalS(double area) const
    {
            return this->diagonal_->computeS(this, area);
    }
    /*! Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *  \param[in] area   area of the i-th tessera to be calculated
     *  \param[in] radius radius of the sphere the tessera belongs to
     */
    virtual double diagonalD(double area, double radius) const
    {
            return this->diagonal_->computeD(this, area, radius);
    }

    virtual void epsilon(double eps) { epsilon_ = eps; }
    virtual double epsilon() const { return epsilon_; }

    friend std::ostream & operator<<(std::ostream & os, UniformDielectric & gf) {
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
        DerivativeTraits distance = sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) +
                          (sp[1] - pp[1]) * (sp[1] - pp[1]) +
                          (sp[2] - pp[2]) * (sp[2] - pp[2]));
        return 1/(epsilon_ * distance);
    }
    virtual std::ostream & printObject(std::ostream & os)
    {
        os << "Green's function type: uniform dielectric" << std::endl;
        os << "Permittivity = " << epsilon_ << std::endl;
        return os;
    }
    double epsilon_;
};

namespace
{
    struct buildUniformDielectric {
        template <typename DerivativeType>
        IGreensFunction * operator()(const greenData & _data) {
            DiagonalIntegrator * integrator = 
		    DiagonalIntegratorFactory::TheDiagonalIntegratorFactory().createDiagonalIntegrator(_data.integratorType);
            return new UniformDielectric<DerivativeType>(_data.epsilon, integrator);
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
