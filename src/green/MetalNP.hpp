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

#ifndef METALNP_HPP
#define METALNP_HPP

#include <complex>
#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

class Element;

#include "GreenData.hpp"
#include "GreensFunction.hpp"

/*! \file MetalNP.hpp
 *  \class MetalNP
 *  \brief Class to describe spherical metal nanoparticles.
 *  \author Stefano Corni, Luca Frediani, Roberto Di Remigio
 *  \date 2011, 2014
 *
 *  This class is a wrapper around the Fortran routines written
 *  by Stefano Corni et al. to take into account the presence
 *  of metal nanoparticles.
 *  References:
 *  http://dx.doi.org/10.1063/1.1342241
 *  http://dx.doi.org/10.1063/1.1507579
 *  http://dx.doi.org/10.1063/1.1558036
 */

template <typename DerivativeTraits,
          typename IntegratorPolicy>
class MetalNP : public GreensFunction<DerivativeTraits, IntegratorPolicy, Metal,
                                      MetalNP<DerivativeTraits, IntegratorPolicy> >
{
private:
    typedef std::complex<double> dcomplex;
public:
    /*! \param[in] e   solvent permittivity
     *  \param[in] eRe real part of the metal nanosphere permittivity
     *  \param[in] eIm imaginary part of the matal nanosphere permittivity
     *  \param[in] o   center of the metal nanosphere
     *  \param[in] r   radius of the metal nanosphere
     */
    MetalNP(double e, double eRe, double eIm, const Eigen::Vector3d & o, double r)
        : GreensFunction<DerivativeTraits, IntegratorPolicy, Metal,
                MetalNP<DerivativeTraits, IntegratorPolicy> >(),
                epsilonSolvent_(e), spherePosition_(o), sphereRadius_(r) { this->profile_ = Metal(eRe, eIm); }
    virtual ~MetalNP() {}
    /*! Returns value of the directional derivative of the
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point, thus assuming that the direction is relative to that point.
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double kernelD(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const override
    {
        return this->epsilonSolvent_ * (this->derivativeProbe(direction, p1, p2));
    }

    /*! Calculates the matrix representation of the S operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd singleLayer(const std::vector<Element> & e) const override
    {
        return this->integrator_.singleLayer(*this, e);
    }
    /*! Calculates the matrix representation of the D operator
     *  \param[in] e list of finite elements
     */
    virtual Eigen::MatrixXd doubleLayer(const std::vector<Element> & e) const override
    {
        return this->integrator_.doubleLayer(*this, e);
    }

    friend std::ostream & operator<<(std::ostream & os, MetalNP & gf) {
        return gf.printObject(os);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    /*! Evaluates the Green's function given a pair of points
     *  \param[in] sp the source point
     *  \param[in] pp the probe point
     */
    virtual DerivativeTraits operator()(DerivativeTraits * sp, DerivativeTraits * pp) const override
    {
        // Just to silence warnings...
        DerivativeTraits distance = sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) +
                          (sp[1] - pp[1]) * (sp[1] - pp[1]) +
                          (sp[2] - pp[2]) * (sp[2] - pp[2]));
        return 1/(this->epsilonSolvent_ * distance);
    }
    double epsilonSolvent_;
    Eigen::Vector3d spherePosition_;
    double sphereRadius_;
    virtual std::ostream & printObject(std::ostream & os)
    {
        os << "Green's function type: metal sphere" << std::endl;
        os << this->profile_ << std::endl;
        os << "Sphere position               = " << spherePosition_ << std::endl;
        os << "Sphere radius                 = " << sphereRadius_ << std::endl;
        os << "Solvent permittivity          = " << epsilonSolvent_;
        return os;
    }
};

#endif // METALNP_HPP
