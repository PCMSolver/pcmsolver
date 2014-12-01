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

#ifndef PLANARINTERFACE_HPP
#define PLANARINTERFACE_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

#include "DerivativeTypes.hpp"
#include "DiagonalIntegrator.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file PlanarInterface.hpp
 *  \class PlanarInterface
 *  \brief Green's functions for a planar interface
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013-2014
 *  \tparam T evaluation strategy for the function and its derivatives
 *
 *  Reference:
 *  http://dx.doi.org/10.1063/1.1643727
 */

template <typename T>
class PlanarInterface : public GreensFunction<T>
{
public:
    PlanarInterface(double eps1, double eps2, const Eigen::Vector3d & pos, double width)
        : GreensFunction<T>(false), eps1_(eps1), eps2_(eps2), pos_(pos), width_(width),
          computed_(false) {}
    PlanarInterface(double eps1, double eps2, const Eigen::Vector3d & pos, double width, DiagonalIntegrator * diag)
        : GreensFunction<T>(false, diag), eps1_(eps1), eps2_(eps2), pos_(pos), width_(width),
          computed_(false) {}
    virtual ~PlanarInterface() {}
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
     *  \param[in] area   area of the i-th tessera to be calculated
     */
    virtual double diagonalS(double area) const;
    /*!
     *  Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *  \param[in] area   area of the i-th tessera to be calculated
     *  \param[in] radius radius of the sphere the tessera belongs to
     */
    virtual double diagonalD(double area, double radius) const;

    virtual double epsilon() const { return eps1_; } // This is just to get it to compile...

    friend std::ostream & operator<<(std::ostream & os, PlanarInterface & gf) {
        return gf.printObject(os);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    /*!
     *  Evaluates the Green's function given a pair of points
     *
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual T operator()(T * source, T * probe) const;
    double eps1_;
    double eps2_;
    Eigen::Vector3d pos_;
    double width_;
    bool computed_;
    virtual std::ostream & printObject(std::ostream & os);
};

namespace
{
    struct buildPlanarInterface {
        template <typename DerivativeType>
        IGreensFunction * operator()(const greenData & _data) {
            // We pass some bogus arguments...
            Eigen::Vector3d orig;
            orig << 0.0, 0.0, 0.0;
            return new PlanarInterface<DerivativeType>(_data.epsilon, 0.0, orig, 0.0, _data.integrator);
        }
    };

    IGreensFunction * createPlanarInterface(const greenData & _data)
    {
        buildPlanarInterface build;
        return for_id<derivative_types>(build, _data, _data.how);
    }
    const std::string PLANARINTERFACE("PLANARINTERFACE");
    const bool registeredPlanarInterface =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            PLANARINTERFACE, createPlanarInterface);
}

#endif // PLANARINTERFACE_HPP
