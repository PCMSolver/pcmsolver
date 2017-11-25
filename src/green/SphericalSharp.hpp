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

#include <cmath>
#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "DerivativeUtils.hpp"
#include "GreensFunction.hpp"
#include "bi_operators/IntegratorForward.hpp"
#include "dielectric_profile/Sharp.hpp"
#include "utils/MathUtils.hpp"
#include "utils/legendre.h"

/*! \file SphericalSharp.hpp
 *  \class SphericalSharp
 *  \brief Green's function for a sharp interface with spherical symmetry
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation
 * of S and D
 */

template <typename DerivativeTraits = AD_directional,
          typename IntegratorPolicy = CollocationIntegrator>
class SphericalSharp __final
    : public GreensFunction<DerivativeTraits,
                            IntegratorPolicy,
                            Sharp,
                            SphericalSharp<DerivativeTraits, IntegratorPolicy> > {
public:
  /*! Constructor for a one-layer interface
   * \param[in] e permittivity of the sphere
   * \param[in] esolv permittivity of the solvent
   * \param[in] r radius of the dielectric sphere
   * \param[in] o center of the sphere
   */
  SphericalSharp(double e, double esolv, double r, const Eigen::Vector3d & o)
      : GreensFunction<DerivativeTraits,
                       IntegratorPolicy,
                       Sharp,
                       SphericalSharp<DerivativeTraits, IntegratorPolicy> >(),
        origin_(o) {
    this->profile_ = Sharp(e, esolv, r);
  }
  /*! Constructor for a one-layer interface
   * \param[in] e permittivity of the sphere
   * \param[in] esolv permittivity of the solvent
   * \param[in] r radius of the dielectric sphere
   * \param[in] o center of the sphere
   */
  SphericalSharp(double e,
                 double esolv,
                 double r,
                 const Eigen::Vector3d & o,
                 double f)
      : GreensFunction<DerivativeTraits,
                       IntegratorPolicy,
                       Sharp,
                       SphericalSharp<DerivativeTraits, IntegratorPolicy> >(f),
        origin_(o) {
    this->profile_ = Sharp(e, esolv, r);
  }
  virtual ~SphericalSharp() {}

  /*! Calculates the matrix representation of the S operator
   *  \param[in] e list of finite elements
   */
  virtual Eigen::MatrixXd singleLayer(const std::vector<Element> & e) const
      __override {
    return this->integrator_.singleLayer(*this, e);
  }
  /*! Calculates the matrix representation of the D operator
   *  \param[in] e list of finite elements
   */
  virtual Eigen::MatrixXd doubleLayer(const std::vector<Element> & e) const
      __override {
    return this->integrator_.doubleLayer(*this, e);
  }

  /*! \brief Returns non-singular part of the Green's function (image potential)
   *  \param[in] source location of the source charge
   *  \param[in] probe location of the probe charge
   */
  double imagePotential(const Eigen::Vector3d & source,
                        const Eigen::Vector3d & probe) const {
    DerivativeTraits sp[3], pp[3];
    sp[0] = source(0);
    sp[1] = source(1);
    sp[2] = source(2);
    pp[0] = probe(0);
    pp[1] = probe(1);
    pp[2] = probe(2);
    return this->imagePotential_impl(sp, pp)[0];
  }
  /*! Returns value of the directional derivative of the
   *  non-singular part (image potential) of the Greens's function for the pair of
   *points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
   *\mathbf{n}_{\mathbf{p}_2}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the probe point, thus assuming that the direction is relative to that point.
   *
   *  \param[in] direction the direction
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   */
  double imagePotentialDerivative(const Eigen::Vector3d & direction,
                                  const Eigen::Vector3d & p1,
                                  const Eigen::Vector3d & p2) const {
    DerivativeTraits t1[3], t2[3];
    t1[0] = p1(0);
    t1[1] = p1(1);
    t1[2] = p1(2);
    t2[0] = p2(0);
    t2[1] = p2(1);
    t2[2] = p2(2);
    t2[0][1] = direction(0);
    t2[1][1] = direction(1);
    t2[2][1] = direction(2);
    return this->imagePotential_impl(t1, t2)[1];
  }

  friend std::ostream & operator<<(std::ostream & os, SphericalSharp & gf) {
    return gf.printObject(os);
  }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See
                                     http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
                                     */
      private :
      /*! Center of the dielectric sphere */
      Eigen::Vector3d origin_;

  /*! Evaluates the Green's function given a pair of points
   *  \param[in] sp the source point
   *  \param[in] pp the probe point
   */
  virtual DerivativeTraits operator()(DerivativeTraits * sp,
                                      DerivativeTraits * pp) const {
    DerivativeTraits distance =
        sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) + (sp[1] - pp[1]) * (sp[1] - pp[1]) +
             (sp[2] - pp[2]) * (sp[2] - pp[2]));
    return (1 / (this->profile_.epsilonSolvent * distance) -
            this->imagePotential_impl(sp, pp));
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
   */
  virtual double kernelD_impl(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const {
    return (this->profile_.epsilonSolvent *
            this->derivativeProbe(direction, p1, p2));
  }
  virtual KernelS exportKernelS_impl() const __override {
    return pcm::bind(&SphericalSharp<DerivativeTraits, IntegratorPolicy>::kernelS,
                     *this,
                     pcm::_1,
                     pcm::_2);
  }
  virtual KernelD exportKernelD_impl() const __override {
    return pcm::bind(&SphericalSharp<DerivativeTraits, IntegratorPolicy>::kernelD,
                     *this,
                     pcm::_1,
                     pcm::_2,
                     pcm::_3);
  }
  DerivativeTraits imagePotential_impl(DerivativeTraits * sp,
                                       DerivativeTraits * pp) const {
    // Data from permittivity profile
    double radius = this->profile_.radius;
    double epsSolv = this->profile_.epsilonSolvent;
    double eps = this->profile_.epsilon;
    DerivativeTraits sp_origin[3], pp_origin[3];
    for (int i = 0; i < 3; ++i) {
      sp_origin[i] = sp[i] - origin_(i);
      pp_origin[i] = pp[i] - origin_(i);
    }
    DerivativeTraits sp_origin_norm = norm(sp_origin);
    DerivativeTraits pp_origin_norm = norm(pp_origin);
    // Angle between source and probe point
    DerivativeTraits cos_gamma =
        dot_product(sp_origin, pp_origin) / (sp_origin_norm * pp_origin_norm);

    DerivativeTraits r_img[3];
    for (int i = 0; i < 3; ++i) {
      r_img[i] = origin_(i) + pow(radius / pp_origin_norm, 2) * pp_origin[i];
    }
    // Distance between sp and r_img
    DerivativeTraits sp_image = distance(sp, r_img);
    DerivativeTraits q_img = radius / pp_origin_norm;
    double factor = (eps - epsSolv) / (eps + epsSolv);
    // DerivativeTraits G_img = DerivativeTraits(0.0);
    DerivativeTraits G_img = factor * (q_img / sp_image - q_img / sp_origin_norm);
    DerivativeTraits f_0 = radius / (sp_origin_norm * pp_origin_norm);
    DerivativeTraits f_L = f_0;
    for (int L = 1; L <= 200; ++L) {
      f_L = f_L * radius * f_0;
      double C_0_L = (eps - epsSolv) * L / ((eps + epsSolv) * L + epsSolv);
      DerivativeTraits pl_x = Legendre::Pn<DerivativeTraits>(L, cos_gamma);
      G_img += f_L * (C_0_L - factor) * pl_x;
    }

    return G_img / epsSolv;
  }
  virtual std::ostream & printObject(std::ostream & os) {
    Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, 0, ", ", "\n", "(", ")");
    os << "Green's function type: spherical sharp" << std::endl;
    os << this->profile_ << std::endl;
    os << "Sphere center        = " << this->origin_.transpose().format(CleanFmt);
    return os;
  }
};
