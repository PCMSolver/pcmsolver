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
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#ifndef IGREENSFUNCTION_HPP
#define IGREENSFUNCTION_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Dense>

class DiagonalIntegrator;

/*! \file IGreensFunction.hpp
 *  \class IGreensFunction
 *  \brief Interface for Green's function classes
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2014
 */

class IGreensFunction
{
public:
    IGreensFunction() : uniform_(false) {}
    IGreensFunction(bool uniform) : uniform_(uniform), diagonal_(NULL) {}
    IGreensFunction(bool uniform, DiagonalIntegrator * diag) : uniform_(uniform), diagonal_(diag) {}
    virtual ~IGreensFunction() {}
    /*!
     *  Returns value of the Greens's function for the pair
     *  of points p1, p2: \f$ G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    virtual double function(const Eigen::Vector3d & p1,
                            const Eigen::Vector3d &p2) const = 0;
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
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const = 0;
    /*!
     *  Returns value of the directional derivative of the
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
                                    const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const = 0;
    /*!
     *  Returns value of the directional derivative of the
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point.
     *
     *  \param[in] normal_p2 the normal vector to p1
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivativeProbe(const Eigen::Vector3d & normal_p2,
                                   const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const = 0;
    /*!
     *  Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the source point.
     *
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    virtual Eigen::Vector3d gradientSource(const Eigen::Vector3d & p1,
                                           const Eigen::Vector3d & p2) const = 0;
    /*!
     *  Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the probe point.
     *
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    virtual Eigen::Vector3d gradientProbe(const Eigen::Vector3d & p1,
                                          const Eigen::Vector3d & p2) const = 0;
    /*!
     *  Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the source point.
     *
     *  \param[in] gradient the gradient
     *  \param[in]       p1 first point
     *  \param[in]       p2 second point
     */
    virtual void gradientSource(Eigen::Vector3d & gradient, const Eigen::Vector3d & p1,
                                const Eigen::Vector3d & p2) const = 0;
    /*!
     *  Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the probe point.
     *
     *  \param[in] gradient the gradient
     *  \param[in]       p1 first point
     *  \param[in]       p2 second point
     */
    virtual void gradientProbe(Eigen::Vector3d & gradient, const Eigen::Vector3d & p1,
                               const Eigen::Vector3d & p2) const = 0;

    /*!
     *  Calculates the diagonal elements of the S operator: \f$ S_{ii} \f$
     *  \param[in] area   area of the i-th tessera to be calculated
     */
    virtual double diagonalS(double area) const = 0;
    /*!
     *  Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *  \param[in] area   area of the i-th tessera to be calculated
     */
    virtual double diagonalD(double area, double radius) const = 0;

    virtual double epsilon() const = 0;
    bool isUniform() const { return uniform_; }

    friend std::ostream & operator<<(std::ostream & os, IGreensFunction & gf) {
        return gf.printObject(os);
    }
protected:
    virtual std::ostream & printObject(std::ostream & os) = 0;
    bool uniform_;
    DiagonalIntegrator * diagonal_;
};

#endif // IGREENSFUNCTION_HPP
