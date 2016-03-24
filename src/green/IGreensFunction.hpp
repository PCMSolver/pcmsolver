/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#ifndef IGREENSFUNCTION_HPP
#define IGREENSFUNCTION_HPP

#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

class Element;

#include "dielectric_profile/ProfileTypes.hpp"

/*! \file IGreensFunction.hpp
 *  \class IGreensFunction
 *  \brief Interface for Green's function classes
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2015
 *
 *  The Non-Virtual Interface (NVI) idiom is used. Notice also that some of the return
 *  types are in some cases "auto", meaning that we will let the compiler deduce
 *  them.
 */

/*! \typedef KernelS
 *  \brief functor handle to the kernelS method
 */
typedef pcm::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> KernelS;

/*! \typedef KernelD
 *  \brief functor handle to the kernelD method
 */
typedef pcm::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &, const Eigen::Vector3d &)> KernelD;

class IGreensFunction
{
  public:
    virtual ~IGreensFunction() {}
    /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e. the value of the
     *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    double kernelS(const Eigen::Vector3d & p1, const Eigen::Vector3d &p2) const {
      return kernelS_impl(p1, p2);
    }
    /*! Returns value of the kernel of the \f$\mathcal{D}\f$ integral operator for the pair of points p1, p2:
     *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this methods with \f$\mathbf{p}_1\f$
     *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} = \mathbf{n}_{\mathbf{p}_1}\f$
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    double kernelD(const Eigen::Vector3d & direction, const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const {
      return kernelD_impl(direction, p1, p2);
    }

    KernelS exportKernelS() const {
      return exportKernelS_impl();
    }
    KernelD exportKernelD() const {
      return exportKernelD_impl();
    }

    /*! Whether the Green's function describes a uniform environment */
    virtual bool uniform() const = 0;
    /*! Returns a dielectric permittivity profile */
    virtual Permittivity permittivity() const = 0;

    /*! Calculates the matrix representation of the S operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd singleLayer(const std::vector<Element> & e) const = 0;
    /*! Calculates the matrix representation of the D operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd doubleLayer(const std::vector<Element> & e) const = 0;

    friend std::ostream & operator<<(std::ostream & os, IGreensFunction & gf) {
      return gf.printObject(os);
    }
  protected:
    /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e. the value of the
     *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    virtual double kernelS_impl(const Eigen::Vector3d & p1, const Eigen::Vector3d &p2) const = 0;
    /*! Returns value of the kernel of the \f$\mathcal{D}\f$ integral operator for the pair of points p1, p2:
     *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this methods with \f$\mathbf{p}_1\f$
     *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} = \mathbf{n}_{\mathbf{p}_1}\f$
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double kernelD_impl(const Eigen::Vector3d & direction,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const = 0;
    virtual KernelS exportKernelS_impl() const = 0;
    virtual KernelD exportKernelD_impl() const = 0;
    virtual std::ostream & printObject(std::ostream & os) = 0;
};

#endif // IGREENSFUNCTION_HPP
