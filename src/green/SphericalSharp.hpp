/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
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

#include <cmath>
#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

/*! \file SphericalSharp.hpp */

namespace pcm {
namespace cavity {
class Element;
} // namespace cavity
namespace dielectric_profile {
struct Sharp;
} // namespace dielectric_profile
} // namespace pcm

#include "DerivativeTypes.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "dielectric_profile/Sharp.hpp"

namespace pcm {
namespace green {
/*! \class SphericalSharp
 *  \brief Green's function for a sharp interface with spherical symmetry
 *  \author Roberto Di Remigio
 *  \date 2018
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 */
template <typename DerivativeTraits = AD_directional>
class SphericalSharp __final
    : public GreensFunction<DerivativeTraits, dielectric_profile::Sharp> {
public:
  /*! Constructor for a one-layer interface
   * \param[in] e permittivity of the sphere
   * \param[in] esolv permittivity of the solvent
   * \param[in] r radius of the dielectric sphere
   * \param[in] o center of the sphere
   * \param[in] l maximum value of angular momentum
   */
  SphericalSharp(double e, double esolv, double r, const Eigen::Vector3d & o, int l);

  virtual double permittivity() const __override __final {
    PCMSOLVER_ERROR("permittivity() only implemented for uniform dielectrics");
  }

  /*! \brief Returns non-singular part of the Green's function (image potential)
   *  \param[in] source location of the source charge
   *  \param[in] probe location of the probe charge
   */
  double imagePotential(const Eigen::Vector3d & source,
                        const Eigen::Vector3d & probe) const;

  /*! Returns value of the directional derivative of the
   *  non-singular part (image potential) of the Greens's function for the pair of
   *  points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1,
   * \mathbf{p}_2)\cdot\mathbf{n}_{\mathbf{p}_2}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the probe point, thus assuming that the direction is relative to that point.
   *
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   */
  double imagePotentialDerivative(const Eigen::Vector3d & direction,
                                  const Eigen::Vector3d & p1,
                                  const Eigen::Vector3d & p2) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  /*! Center of the dielectric sphere */
  Eigen::Vector3d origin_;

  /*! Maximum angular momentum in the final summation over Legendre polynomials */
  int maxLGreen_;

  /*! Evaluates the Green's function given a pair of points
   *  \param[in] sp the source point
   *  \param[in] pp the probe point
   */
  virtual DerivativeTraits operator()(DerivativeTraits * sp,
                                      DerivativeTraits * pp) const __override;

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
                              const Eigen::Vector3d & p2) const __override;

  virtual KernelS exportKernelS_impl() const __override;
  virtual KernelD exportKernelD_impl() const __override;
  virtual DerivativeProbe exportDerivativeProbe_impl() const __override;

  virtual double singleLayer_impl(const Element & e, double factor) const __override;
  virtual double doubleLayer_impl(const Element & e, double factor) const __override;

  DerivativeTraits imagePotential_impl(DerivativeTraits * sp,
                                       DerivativeTraits * pp) const;

  virtual std::ostream & printObject(std::ostream & os) __override;
};

template <typename DerivativeTraits>
IGreensFunction * createSphericalSharp(const GreenData & data) {
  return new SphericalSharp<DerivativeTraits>(
      data.epsilon1, data.epsilon2, data.center, data.origin, data.maxL);
}
} // namespace green
} // namespace pcm
