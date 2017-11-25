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

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "IGreensFunction.hpp"
#include "utils/Stencils.hpp"

/*! \file GreensFunction.hpp */

namespace pcm {
namespace green {
template <typename DerivativeTraits, typename ProfilePolicy>

/*! \class GreensFunction
 *  \brief Templated interface for Green's functions
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2016
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam ProfilePolicy    dielectric profile type
 */
class GreensFunction : public IGreensFunction {
public:
  GreensFunction() : delta_(1.0e-04) {}
  virtual ~GreensFunction() {}
  /*! @{ Methods to sample the Green's function directional derivatives */
  /*! Returns value of the directional derivative of the
   *  Greens's function for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
   * \mathbf{n}_{\mathbf{p}_1}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the source point.
   *  \param[in] normal_p1 the normal vector to p1
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   */
  virtual double derivativeSource(const Eigen::Vector3d & normal_p1,
                                  const Eigen::Vector3d & p1,
                                  const Eigen::Vector3d & p2) const {
    DerivativeTraits t1[3], t2[3];
    t1[0] = p1(0);
    t1[1] = p1(1);
    t1[2] = p1(2);
    t1[0][1] = normal_p1(0);
    t1[1][1] = normal_p1(1);
    t1[2][1] = normal_p1(2);
    t2[0] = p2(0);
    t2[1] = p2(1);
    t2[2] = p2(2);
    return this->operator()(t1, t2)[1];
  }
  /*! Returns value of the directional derivative of the
   *  Greens's function for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
   * \mathbf{n}_{\mathbf{p}_2}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the probe point.
   *  \param[in] normal_p2 the normal vector to p2
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   */
  virtual double derivativeProbe(const Eigen::Vector3d & normal_p2,
                                 const Eigen::Vector3d & p1,
                                 const Eigen::Vector3d & p2) const {
    DerivativeTraits t1[3], t2[3];
    t1[0] = p1(0);
    t1[1] = p1(1);
    t1[2] = p1(2);
    t2[0] = p2(0);
    t2[1] = p2(1);
    t2[2] = p2(2);
    t2[0][1] = normal_p2(0);
    t2[1][1] = normal_p2(1);
    t2[2][1] = normal_p2(2);
    return this->operator()(t1, t2)[1];
  }
  /*! @} */

  /*! @{ Methods to sample the Green's function gradients */
  /*! Returns full gradient of Greens's function for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
   *  Notice that this method returns the gradient with respect to the source point.
   *  \param[in] p1 first point
   *  \param[in] p2 second point
   */
  Eigen::Vector3d gradientSource(const Eigen::Vector3d & p1,
                                 const Eigen::Vector3d & p2) const {
    return (Eigen::Vector3d() << derivativeSource(Eigen::Vector3d::UnitX(), p1, p2),
            derivativeSource(Eigen::Vector3d::UnitY(), p1, p2),
            derivativeSource(Eigen::Vector3d::UnitZ(), p1, p2))
        .finished();
  }
  /*! Returns full gradient of Greens's function for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
   *  Notice that this method returns the gradient with respect to the probe point.
   *  \param[in] p1 first point
   *  \param[in] p2 second point
   */
  Eigen::Vector3d gradientProbe(const Eigen::Vector3d & p1,
                                const Eigen::Vector3d & p2) const {
    return (Eigen::Vector3d() << derivativeProbe(Eigen::Vector3d::UnitX(), p1, p2),
            derivativeProbe(Eigen::Vector3d::UnitY(), p1, p2),
            derivativeProbe(Eigen::Vector3d::UnitZ(), p1, p2))
        .finished();
  }
  /*! @}*/

  /*! Whether the Green's function describes a uniform environment */
  virtual bool uniform() const __final __override {
    return dielectric_profile::uniform(this->profile_);
  }
  /*! Returns a dielectric permittivity profile */
  virtual Permittivity permittivity() const __final __override {
    return this->profile_;
  }

  friend std::ostream & operator<<(std::ostream & os, GreensFunction & gf) {
    return gf.printObject(os);
  }

protected:
  /*! Evaluates the Green's function given a pair of points
   *  \param[in] source the source point
   *  \param[in]  probe the probe point
   */
  virtual DerivativeTraits operator()(DerivativeTraits * source,
                                      DerivativeTraits * probe) const = 0;
  /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e.
   *  the value of the Greens's function for the pair of points p1, p2:
   *  \f$ G(\mathbf{p}_1, \mathbf{p}_2)\f$
   *  \param[in] p1 first point
   *  \param[in] p2 second point
   *  \note Relies on the implementation of operator() in the subclasses and
   *  that is all subclasses need to implement.
   *  Thus this method is marked __final.
   */
  virtual double kernelS_impl(const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const __final __override {
    DerivativeTraits sp[3], pp[3];
    sp[0] = p1(0);
    sp[1] = p1(1);
    sp[2] = p1(2);
    pp[0] = p2(0);
    pp[1] = p2(1);
    pp[2] = p2(2);
    return this->operator()(sp, pp)[0];
  }
  virtual std::ostream & printObject(std::ostream & os) __override {
    os << "Green's Function" << std::endl;
    return os;
  }
  /// Step for numerical differentiation
  double delta_;
  /// Permittivity profile
  ProfilePolicy profile_;
};

template <typename ProfilePolicy>
class GreensFunction<Stencil, ProfilePolicy> : public IGreensFunction {
public:
  GreensFunction() : delta_(1.0e-04) {}
  virtual ~GreensFunction() {}
  /*! @{ Methods to sample the Green's function directional derivatives */
  /*! Returns value of the directional derivative of the
   *  Greens's function for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
   * \mathbf{n}_{\mathbf{p}_1}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the source point.
   *  \param[in] normal_p1 the normal vector to p1
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   */
  virtual double derivativeSource(const Eigen::Vector3d & normal_p1,
                                  const Eigen::Vector3d & p1,
                                  const Eigen::Vector3d & p2) const {
    return threePointStencil(
        pcm::bind(&GreensFunction<Stencil, ProfilePolicy>::kernelS,
                  this,
                  pcm::_1,
                  pcm::_2),
        p1,
        p2,
        normal_p1,
        this->delta_);
  }
  /*! Returns value of the directional derivative of the
   *  Greens's function for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot
   * \mathbf{n}_{\mathbf{p}_2}\f$
   *  Notice that this method returns the directional derivative with respect
   *  to the probe point.
   *  \param[in] normal_p2 the normal vector to p2
   *  \param[in]        p1 first point
   *  \param[in]        p2 second point
   */
  virtual double derivativeProbe(const Eigen::Vector3d & normal_p2,
                                 const Eigen::Vector3d & p1,
                                 const Eigen::Vector3d & p2) const __final {
    return threePointStencil(
        pcm::bind(&GreensFunction<Stencil, ProfilePolicy>::kernelS,
                  this,
                  pcm::_1,
                  pcm::_2),
        p2,
        p1,
        normal_p2,
        this->delta_);
  }
  /*! @}*/

  /*! @{ Methods to sample the Green's functions gradients */
  /*! Returns full gradient of Greens's function for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
   *  Notice that this method returns the gradient with respect to the source point.
   *  \param[in] p1 first point
   *  \param[in] p2 second point
   */
  Eigen::Vector3d gradientSource(const Eigen::Vector3d & p1,
                                 const Eigen::Vector3d & p2) const {
    return (Eigen::Vector3d() << derivativeSource(Eigen::Vector3d::UnitX(), p1, p2),
            derivativeSource(Eigen::Vector3d::UnitY(), p1, p2),
            derivativeSource(Eigen::Vector3d::UnitZ(), p1, p2))
        .finished();
  }
  /*! Returns full gradient of Greens's function for the pair of points p1, p2:
   *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
   *  Notice that this method returns the gradient with respect to the probe point.
   *  \param[in] p1 first point
   *  \param[in] p2 second point
   */
  Eigen::Vector3d gradientProbe(const Eigen::Vector3d & p1,
                                const Eigen::Vector3d & p2) const {
    return (Eigen::Vector3d() << derivativeProbe(Eigen::Vector3d::UnitX(), p1, p2),
            derivativeProbe(Eigen::Vector3d::UnitY(), p1, p2),
            derivativeProbe(Eigen::Vector3d::UnitZ(), p1, p2))
        .finished();
  }
  /*! @}*/

  /*! Whether the Green's function describes a uniform environment */
  virtual bool uniform() const __final __override {
    return dielectric_profile::uniform(this->profile_);
  }
  /*! Returns a dielectric permittivity profile */
  virtual Permittivity permittivity() const __final __override {
    return this->profile_;
  }

  friend std::ostream & operator<<(std::ostream & os, GreensFunction & gf) {
    return gf.printObject(os);
  }

protected:
  /*! Evaluates the Green's function given a pair of points
   *  \param[in] source the source point
   *  \param[in]  probe the probe point
   */
  virtual Stencil operator()(Stencil * source, Stencil * probe) const = 0;
  /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e.
   * the value of the
   *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1,
   * \mathbf{p}_2)\f$
   *  \param[in] p1 first point
   *  \param[in] p2 second point
   *  \note Relies on the implementation of operator() in the subclasses and that is
   * all subclasses
   *  need to implement. Thus this method is marked __final.
   */
  virtual double kernelS_impl(const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const __final __override {
    return this->operator()(const_cast<Stencil *>(p1.data()),
                            const_cast<Stencil *>(p2.data()));
  }
  virtual std::ostream & printObject(std::ostream & os) __override {
    os << "Green's Function" << std::endl;
    return os;
  }
  /// Step for numerical differentiation
  double delta_;
  /// Permittivity profile
  ProfilePolicy profile_;
};

namespace detail {
/*! Approximate collocation formula for the diagonal of:
 * \f[
 *  (\mathcal{S}_\mathrm{i}f)(\mathbf{s}) = \int_\Gamma
 * \mathop{}\!\mathrm{d}\mathbf{s}^\prime
 *  \frac{f(\mathbf{s}^\prime)}{|\mathbf{s} - \mathbf{s}^\prime|}
 *  \f]
 *  The
 *  \param[in] area   area of the finite element, \f$ a_i \f$
 *  \param[in] factor scaling factor for diagonal elements, \f$ k \f$
 *  \return The approximate value of \f$S_{ii, \mathrm{i}} = k\sqrt{\frac{4\pi}{a_i}}
 * \f$
 */
inline double diagonalSi(double area, double factor) {
  return (factor * std::sqrt(4 * M_PI / area));
}

/*! Approximate collocation formula for the diagonal of:
 * \f[
 *  (\mathcal{D}_\mathrm{i}f)(\mathbf{s}) = \int_\Gamma
 * \mathop{}\!\mathrm{d}\mathbf{s}^\prime
 *  \frac{(\mathbf{s} - \mathbf{s}^\prime)\cdot
 * \mathbf{n}_{\mathbf{s}^\prime}}{|\mathbf{s} -
 * \mathbf{s}^\prime|^3}f(\mathbf{s}^\prime)
 *  \f]
 *  The
 *  \param[in] area   area of the finite element, \f$ a_i \f$
 *  \param[in] radius radius of the finite element, \f$ R_i \f$
 *  \param[in] factor scaling factor for diagonal elements, \f$ k \f$
 *  \return The approximate value of \f$D_{ii, \mathrm{i}} =
 * -k\sqrt{\frac{\pi}{a_i}}\frac{1}{R_i} \f$
 */
inline double diagonalDi(double area, double radius, double factor) {
  return (-factor * std::sqrt(M_PI / area) * 1.0 / radius);
}
} // namespace detail
} // namespace green
} // namespace pcm
