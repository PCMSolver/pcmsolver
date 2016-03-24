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

#ifndef IONICLIQUID_HPP
#define IONICLIQUID_HPP

#include <cmath>
#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

class Element;

#include "DerivativeTypes.hpp"
#include "DerivativeUtils.hpp"
#include "bi_operators/IntegratorForward.hpp"
#include "GreensFunction.hpp"
#include "dielectric_profile/Yukawa.hpp"

/*! \file IonicLiquid.hpp
 *  \class IonicLiquid
 *  \brief Green's functions for ionic liquid, described by the linearized Poisson-Boltzmann equation.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013-2015
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation of S and D
 */

template <typename DerivativeTraits = AD_directional,
          typename IntegratorPolicy = CollocationIntegrator>
class IonicLiquid __final : public GreensFunction<DerivativeTraits, IntegratorPolicy, Yukawa,
                                     IonicLiquid<DerivativeTraits, IntegratorPolicy> >
{
public:
    IonicLiquid(double eps, double k) : GreensFunction<DerivativeTraits, IntegratorPolicy, Yukawa,
                                                  IonicLiquid<DerivativeTraits, IntegratorPolicy> >() { this->profile_ = Yukawa(eps, k); }
    IonicLiquid(double eps, double k, double f) : GreensFunction<DerivativeTraits, IntegratorPolicy, Yukawa,
                                                  IonicLiquid<DerivativeTraits, IntegratorPolicy> >(f) { this->profile_ = Yukawa(eps, k); }
    virtual ~IonicLiquid() {}

    /*! Calculates the matrix representation of the S operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd singleLayer(const std::vector<Element> & e) const __override
    {
        return this->integrator_.singleLayer(*this, e);
    }
    /*! Calculates the matrix representation of the D operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd doubleLayer(const std::vector<Element> & e) const __override
    {
        return this->integrator_.doubleLayer(*this, e);
    }

    friend std::ostream & operator<<(std::ostream & os, IonicLiquid & gf) {
        return gf.printObject(os);
    }
private:
    /*!
     *  Evaluates the Green's function given a pair of points
     *
     *  \param[in] sp the source point
     *  \param[in] pp the probe point
     */
    virtual DerivativeTraits operator()(DerivativeTraits * sp, DerivativeTraits * pp) const __override
    {
        double eps = this->profile_.epsilon;
	    double k = this->profile_.kappa;
        return (exp(-k * distance(sp, pp)) / (eps * distance(sp, pp)));
    }
    /*! Returns value of the directional derivative of the
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point, thus assuming that the direction is relative to that point.
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double kernelD_impl(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const __override
    {
        return this->profile_.epsilon * (this->derivativeProbe(direction, p1, p2));
    }
    virtual KernelS exportKernelS_impl() const __override {
      return pcm::bind(&IonicLiquid<DerivativeTraits, IntegratorPolicy>::kernelS, *this, pcm::_1, pcm::_2);
    }
    virtual KernelD exportKernelD_impl() const __override {
      return pcm::bind(&IonicLiquid<DerivativeTraits, IntegratorPolicy>::kernelD, *this, pcm::_1, pcm::_2, pcm::_3);
    }
    virtual std::ostream & printObject(std::ostream & os) __override
    {
        os << "Green's function type: ionic liquid" << std::endl;
        os << this->profile_;
        return os;
    }
};

#endif // IONICLIQUID_HPP
