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

#ifndef UNIFORMDIELECTRIC_HPP
#define UNIFORMDIELECTRIC_HPP

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
#include "dielectric_profile/Uniform.hpp"

/*! \file UniformDielectric.hpp
 *  \class UniformDielectric
 *  \brief Green's function for uniform dielectric.
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2015
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation of S and D
 */

template <typename DerivativeTraits = AD_directional,
          typename IntegratorPolicy = CollocationIntegrator>
class UniformDielectric __final : public GreensFunction<DerivativeTraits, IntegratorPolicy, Uniform,
                                     UniformDielectric<DerivativeTraits, IntegratorPolicy> >
{
public:
    UniformDielectric(double eps) : GreensFunction<DerivativeTraits, IntegratorPolicy, Uniform,
                              UniformDielectric<DerivativeTraits, IntegratorPolicy> >() { this->profile_ = Uniform(eps); }
    UniformDielectric(double eps, double f) : GreensFunction<DerivativeTraits, IntegratorPolicy, Uniform,
                              UniformDielectric<DerivativeTraits, IntegratorPolicy> >(f) { this->profile_ = Uniform(eps); }
    virtual ~UniformDielectric() {}

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

    double epsilon() const { return this->profile_.epsilon; }

    friend std::ostream & operator<<(std::ostream & os, UniformDielectric & gf) {
        return gf.printObject(os);
    }
private:
    /*! Evaluates the Green's function given a pair of points
     *
     *  \param[in] sp the source point
     *  \param[in] pp the probe point
     */
    virtual DerivativeTraits operator()(DerivativeTraits * sp, DerivativeTraits * pp) const __override
    {
        return 1/(this->profile_.epsilon * distance(sp, pp));
    }
    /*! Returns value of the kernel of the \f$\mathcal{D}\f$ integral operator
     *  for the pair of points p1, p2:
     *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this methods with \f$\mathbf{p}_1\f$
     *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} = \mathbf{n}_{\mathbf{p}_1}\f$
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
      return pcm::bind(&UniformDielectric<DerivativeTraits, IntegratorPolicy>::kernelS, *this, pcm::_1, pcm::_2);
    }
    virtual KernelD exportKernelD_impl() const __override {
      return pcm::bind(&UniformDielectric<DerivativeTraits, IntegratorPolicy>::kernelD, *this, pcm::_1, pcm::_2, pcm::_3);
    }
    virtual std::ostream & printObject(std::ostream & os) __override
    {
        os << "Green's function type: uniform dielectric" << std::endl;
        os << this->profile_;
        return os;
    }
};

#endif // UNIFORMDIELECTRIC_HPP
