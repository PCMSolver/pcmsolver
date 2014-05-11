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

#ifndef METALNP_HPP
#define METALNP_HPP

#include <complex>
#include <iosfwd>
#include <string>

#include "Config.hpp"

#include "EigenPimpl.hpp"

class Element;

#include "DiagonalIntegrator.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

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

class MetalNP : public GreensFunction<double>
{
private:
    typedef std::complex<double> dcomplex;
public:
    MetalNP(double eps, double epsRe, double epsIm,
                const Eigen::Vector3d & pos, double radius)
        : GreensFunction<double>(false), epsSolvent_(eps), epsMetal_(dcomplex(epsRe, epsIm)),
          sphPosition_(pos), sphRadius_(radius) {}
    MetalNP(double eps, double epsRe, double epsIm,
                const Eigen::Vector3d & pos, double radius, DiagonalIntegrator * diag)
        : GreensFunction<double>(false, diag), epsSolvent_(eps), epsMetal_(dcomplex(epsRe, epsIm)),
          sphPosition_(pos), sphRadius_(radius) {}
    virtual ~MetalNP() {}
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
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;

    virtual double diagonalS(const Element & e) const {
	    return 1.0;
    }
    virtual double diagonalD(const Element & e) const {
	    return 1.0;
    }

    virtual double epsilon() const { return epsSolvent_; } // This is just to get it to compile...

    friend std::ostream & operator<<(std::ostream & os, MetalNP & gf) {
        return gf.printObject(os);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    /*!
     *  Evaluates the Green's function given a pair of points
     *
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual double operator()(double * source, double * probe) const;
    double epsSolvent_;
    dcomplex epsMetal_;
    Eigen::Vector3d sphPosition_;
    double sphRadius_;
    virtual std::ostream & printObject(std::ostream & os);
};

namespace
{
    // The build functor and use of for_id are not necessary as MetalNP
    // inherits from a GreensFunction<double>
    IGreensFunction * createMetalNP(const greenData & _data)
    {
	// The NP center is in a std::vector<double> but we need an Eigen::Vector3d
	// We are currently assuming that there is only one spherical metal NP
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	for (int i = 0; i < 3; ++i) {
		center(i) = _data.NPspheres[i];
	}
        return new MetalNP(_data.epsilon, _data.epsilonReal, _data.epsilonImaginary, center, _data.NPradii, _data.integrator);
    }
    const std::string METALNP("MetalNP");
    const bool registeredMetalNP =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            METALNP, createMetalNP);
}
#endif // METALNP_HPP
