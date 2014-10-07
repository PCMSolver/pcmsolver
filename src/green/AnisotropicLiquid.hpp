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

#ifndef ANISOTROPICLIQUID_HPP
#define ANISOTROPICLIQUID_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>
#include "taylor.hpp"

class DiagonalIntegrator;

#include "DerivativeTypes.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file AnisotropicLiquid.hpp
 *  \class AnisotropicLiquid
 *  \brief Green's functions for anisotropic liquid, described by a tensorial permittivity
 *  \author Roberto Di Remigio
 *  \date 2014
 *  \tparam T evaluation strategy for the function and its derivatives
 */

template <typename T>
class AnisotropicLiquid : public GreensFunction<T>
{
public:
    /*!
     * \param[in] eigen_eps eigenvalues of the permittivity tensors
     * \param[in] euler_ang Euler angles in degrees
     */
    AnisotropicLiquid(const Eigen::Vector3d & eigen_eps, const Eigen::Vector3d & euler_ang) : 
	    GreensFunction<T>(false), epsilonLab_(eigen_eps), eulerAngles_(euler_ang) { this->build(); }
    /*!
     * \param[in] eigen_eps eigenvalues of the permittivity tensors
     * \param[in] euler_ang Euler angles in degrees
     */
    AnisotropicLiquid(const Eigen::Vector3d & eigen_eps, const Eigen::Vector3d & euler_ang, DiagonalIntegrator * diag) : 
	    GreensFunction<T>(false, diag), epsilonLab_(eigen_eps), eulerAngles_(euler_ang) { this->build(); }
    virtual ~AnisotropicLiquid() {}
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

    virtual double epsilon() const { return epsilonLab_(0); }

    friend std::ostream & operator<<(std::ostream & os, AnisotropicLiquid & gf) {
        return gf.printObject(os);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    /*! Initializes some internals
     */
    void build();	    
    /*!
     *  Evaluates the Green's function given a pair of points
     *
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual T operator()(T * source, T * probe) const;
    /// Diagonal of the permittivity tensor in the lab-fixed frame
    Eigen::Vector3d epsilonLab_;
    /// Euler angles (in degrees) relating molecule-fixed and lab-fixed frames
    Eigen::Vector3d eulerAngles_;
    /// Permittivity tensor in molecule-fixed frame
    Eigen::Matrix3d epsilon_;
    /// Inverse of the permittivity tensor in molecule-fixed frame
    Eigen::Matrix3d epsilonInv_;
    /// molecule-fixed to lab-fixed frames rotation matrix
    Eigen::Matrix3d R_;
    /// Determinant of the permittivity tensor
    double detEps_;
    virtual std::ostream & printObject(std::ostream & os);
};

namespace
{
    struct buildAnisotropicLiquid {
        template <typename DerivativeType>
        IGreensFunction * operator()(const greenData & _data) {
            return new AnisotropicLiquid<DerivativeType>(_data.epsilonTensor, _data.eulerAngles, _data.integrator);
        }
    };

    IGreensFunction * createAnisotropicLiquid(const greenData & _data)
    {
        buildAnisotropicLiquid build;
        return for_id<derivative_types>(build, _data, _data.how);
    }
    const std::string ANISOTROPICLIQUID("AnisotropicLiquid");
    const bool registeredAnisotropicLiquid =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            ANISOTROPICLIQUID, createAnisotropicLiquid);
}

#endif // ANISOTROPICLIQUID_HPP
