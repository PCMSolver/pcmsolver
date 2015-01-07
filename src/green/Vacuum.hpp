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

#ifndef VACUUM_HPP
#define VACUUM_HPP

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

/*! \file Vacuum.hpp
 *  \class Vacuum
 *  \brief Green's function for vacuum.
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2014
 *  \tparam T evaluation strategy for the function and its derivatives
 */

template <typename T>
class Vacuum : public GreensFunction<T>
{
public:
    Vacuum() : GreensFunction<T>(true), epsilon_(1.0) {}
    Vacuum(DiagonalIntegrator * diag) : GreensFunction<T>(true, diag), epsilon_(1.0) {}
    virtual ~Vacuum() {}
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
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    
    /*!
     *  Calculates the diagonal elements of the S operator: \f$ S_{ii} \f$
     *  \param[in] area   area of the i-th tessera
     */
    virtual double diagonalS(double area) const;
    /*!
     *  Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *  \param[in] area   area of the i-th tessera
     *  \param[in] radius radius of the sphere the tessera belongs to
     */
    virtual double diagonalD(double area, double radius) const;

    virtual double epsilon() const { return 1.0; }

    friend std::ostream & operator<<(std::ostream & os, Vacuum & gf) {
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
    struct buildVacuum {
        template <typename DerivativeType>
        IGreensFunction * operator()(const greenData & _data) {
            return new Vacuum<DerivativeType>(_data.integrator);
        }
    };

    IGreensFunction * createVacuum(const greenData & _data)
    {
        buildVacuum build;
        return for_id<derivative_types>(build, _data, _data.how);
    }
    const std::string VACUUM("VACUUM");
    const bool registeredVacuum =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            VACUUM, createVacuum);
}

#endif // VACUUM_HPP
