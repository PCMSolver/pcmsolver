/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#pragma once

#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

/*! \file IGreensFunction.hpp */

namespace pcm {
namespace cavity {
class Element;
} // namespace green
} // namespace pcm

#include "dielectric_profile/ProfileTypes.hpp"

namespace pcm {
using cavity::Element;
using dielectric_profile::Permittivity;
/*! \typedef KernelS
 *  \brief functor handle to the kernelS method
 */
typedef function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> KernelS;

/*! \typedef KernelD
 *  \brief functor handle to the kernelD method
 */
typedef pcm::function<double(const Eigen::Vector3d &,
                             const Eigen::Vector3d &,
                             const Eigen::Vector3d &)>
    KernelD;

/*! \typedef DerivativeProbe
 *  \brief functor handle to the derivativeProbe method
 *  \note This is the directional derivative wrt the probe point
 */
typedef pcm::function<double(const Eigen::Vector3d &,
                             const Eigen::Vector3d &,
                             const Eigen::Vector3d &)>
    DerivativeProbe;

/*! \class IGreensFunction
 *  \brief Interface for Green's function classes
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2016
 *
 *  We **define** as _Green's function_ a function:
 *  \f[
 *      G(\mathbf{r}, \mathbf{r}^\prime) : \mathbb{R}^6 \rightarrow \mathbb{R}
 *  \f]
 *  Green's functions and their directional derivatives appear as kernels of
 *  the \f$\mathcal{S}\f$ and \f$\mathcal{D}\f$ integral operators.
 *  Forming the matrix representation of these operators requires performing
 *  integrations over surface finite elements.
 *  Since these Green's functions present a Coulombic divergence, the diagonal
 *  elements of the operators will diverge unless appropriately formulated.
 *  This is possible, but requires **explicit** access to the _subtype_
 *  of this abstract base object.
 *  This justifies the need for the singleLayer and doubleLayer functions.
 *  The code uses the Non-Virtual Interface (NVI) idiom.
 */
class IGreensFunction {
public:
  virtual ~IGreensFunction() {}

  /*! @{ Methods to sample the Green's function and its probe point directional
   * derivative */
  /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e.
   * the value of the
   *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1,
   * \mathbf{p}_2)\f$
   *  \param[in] p1 first point
   *  \param[in] p2 second point
   *  \note This is the Non-Virtual Interface (NVI)
   */
  double kernelS(const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const {
    return kernelS_impl(p1, p2);
  }
  /*! Returns value of the kernel of the \f$\mathcal{D}\f$ integral operator for the
   * pair of points p1, p2:
   *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1,
   * \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
   *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this
   * methods with \f$\mathbf{p}_1\f$
   *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} =
   * \mathbf{n}_{\mathbf{p}_1}\f$
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   *  \note This is the Non-Virtual Interface (NVI)
   */
  double kernelD(const Eigen::Vector3d & direction,
                 const Eigen::Vector3d & p1,
                 const Eigen::Vector3d & p2) const {
    return kernelD_impl(direction, p1, p2);
  }
  /*! @}*/

  KernelS exportKernelS() const { return exportKernelS_impl(); }
  KernelD exportKernelD() const { return exportKernelD_impl(); }
  DerivativeProbe exportDerivativeProbe() const {
    return exportDerivativeProbe_impl();
  }

  /*! Whether the Green's function describes a uniform environment */
  virtual bool uniform() const = 0;
  /*! Returns a dielectric permittivity profile */
  virtual Permittivity permittivity() const = 0;

  /*! @{ Methods to compute the diagonal of the matrix representation of the S and D
   *    operators by approximate collocation. */
  /*! Calculates an element on the diagonal of the matrix representation of the
   * S operator using an approximate collocation formula.
   *  \param[in] e finite element on the cavity
   *  \param[in] factor the scaling factor for the diagonal elements
   *  \note This is the Non-Virtual Interface (NVI)
   */
  double singleLayer(const Element & e, double factor) const {
    return singleLayer_impl(e, factor);
  }
  /*! Calculates an element of the diagonal of the matrix representation of the D
   * operator
   * using an approximate collocation formula.
   *  \param[in] e finite element on the cavity
   *  \param[in] factor the scaling factor for the diagonal elements
   *  \note This is the Non-Virtual Interface (NVI)
   */
  double doubleLayer(const Element & e, double factor) const {
    return doubleLayer_impl(e, factor);
  }
  /*! @}*/

  friend std::ostream & operator<<(std::ostream & os, IGreensFunction & gf) {
    return gf.printObject(os);
  }

protected:
  /*! @{ Methods to sample the Green's function and its probe point directional
   * derivative */
  /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e.
   * the value of the
   *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1,
   * \mathbf{p}_2)\f$
   *  \param[in] p1 first point
   *  \param[in] p2 second point
   */
  virtual double kernelS_impl(const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const = 0;
  /*! Returns value of the kernel of the \f$\mathcal{D}\f$ integral operator for the
   * pair of points p1, p2:
   *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1,
   * \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
   *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this
   * methods with \f$\mathbf{p}_1\f$
   *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} =
   * \mathbf{n}_{\mathbf{p}_1}\f$
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   */
  virtual double kernelD_impl(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const = 0;
  /*! @}*/

  virtual KernelS exportKernelS_impl() const = 0;
  virtual KernelD exportKernelD_impl() const = 0;
  virtual DerivativeProbe exportDerivativeProbe_impl() const = 0;

  /*! @{ Methods to compute the diagonal of the matrix representation of the S and D
   *    operators by approximate collocation. */
  /*! Calculates an element on the diagonal of the matrix representation of the
   * S operator using an approximate collocation formula.
   *  \param[in] e finite element on the cavity
   *  \param[in] factor the scaling factor for the diagonal elements
   */
  virtual double singleLayer_impl(const Element & e, double factor) const = 0;
  /*! Calculates an element of the diagonal of the matrix representation of the D
   * operator
   * using an approximate collocation formula.
   *  \param[in] e finite element on the cavity
   *  \param[in] factor the scaling factor for the diagonal elements
   */
  virtual double doubleLayer_impl(const Element & e, double factor) const = 0;
  /*! @}*/

  virtual std::ostream & printObject(std::ostream & os) = 0;
};
} // namespace pcm
