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

#include <cmath>
#include <functional>
#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Dense>

class Element;

#include "DerivativeTypes.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file AnisotropicLiquid.hpp
 *  \class AnisotropicLiquid
 *  \brief Green's functions for anisotropic liquid, described by a tensorial permittivity
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam IntegratorPolicy policy for the calculation of diagonal elements
 */

template <typename DerivativeTraits,
          typename IntegratorPolicy>
class AnisotropicLiquid : public GreensFunction<DerivativeTraits, IntegratorPolicy, Anisotropic,
                                     AnisotropicLiquid<DerivativeTraits, IntegratorPolicy> >
{
public:
    /*! \param[in] eigen_eps eigenvalues of the permittivity tensors
     *  \param[in] euler_ang Euler angles in degrees
     */
    AnisotropicLiquid(const Eigen::Vector3d & eigen_eps, const Eigen::Vector3d & euler_ang) :
        GreensFunction<DerivativeTraits, IntegratorPolicy, Anisotropic,
                  AnisotropicLiquid<DerivativeTraits, IntegratorPolicy> >() { this->profile_ = Anisotropic(eigen_eps, euler_ang); }
    virtual ~AnisotropicLiquid() {}
    /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e. the value of the
     *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    virtual double function(const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
    {
        DerivativeTraits sp[3], pp[3], res;
        sp[0] = p1(0); sp[1] = p1(1); sp[2] = p1(2);
        pp[0] = p2(0); pp[1] = p2(1); pp[2] = p2(2);
        res = this->operator()(sp, pp);
        return res[0];
    }
    /*! Returns value of the kernel for the calculation of the \f$\mathcal{D}\f$ integral operator
     *  for the pair of points p1, p2:
     *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this methods with \f$\mathbf{p}_1\f$
     *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} = \mathbf{n}_{\mathbf{p}_1}\f$
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivative(const Eigen::Vector3d & /* direction */,
                              const Eigen::Vector3d & /* p1 */, const Eigen::Vector3d & /* p2 */) const
    {
        // Need the full gradient to get the kernel of D and D^\dagger
        /*
        Eigen::Vector3d scratch = this->profile_.epsilon() * (this->gradientProbe(p1, p2));
        return scratch.dot(direction);
        */
        return 0.0;
    }
    /*! Returns value of the directional derivative of the
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_1}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the source point.
     *  \param[in] normal_p1 the normal vector to p1
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivativeSource(const Eigen::Vector3d & normal_p1,
                            const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
    {
        DerivativeTraits t1[3], t2[3], der;
        t1[0] = p1(0); t1[1] = p1(1); t1[2] = p1(2);
        t1[0][1] = normal_p1(0); t1[1][1] = normal_p1(1); t1[2][1] = normal_p1(2);
        t2[0] = p2(0); t2[1] = p2(1); t2[2] = p2(2);
        der = this->operator()(t1, t2);
        return der[1];
    }
    /*! Returns value of the directional derivative of the
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point.
     *  \param[in] normal_p2 the normal vector to p2
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivativeProbe(const Eigen::Vector3d & normal_p2,
                                   const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
    {
        DerivativeTraits t1[3], t2[3], der;
        t1[0] = p1(0); t1[1] = p1(1); t1[2] = p1(2);
        t2[0] = p2(0); t2[1] = p2(1); t2[2] = p2(2);
        t2[0][1] = normal_p2(0); t2[1][1] = normal_p2(1); t2[2][1] = normal_p2(2);
        der = this->operator()(t1, t2);
        return der[1];
    }

    /*!
     *  Calculates the diagonal elements of the S operator: \f$ S_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalS(const Element & e) const
    {
        return this->diagonal_.computeS(*this, e);
    }
    /*!
     *  Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalD(const Element & e) const
    {
        return this->diagonal_.computeD(*this, e);
    }

    virtual double epsilon() const { return 0.0; }

    friend std::ostream & operator<<(std::ostream & os, AnisotropicLiquid & gf) {
        return gf.printObject(os);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    /*!
     *  Evaluates the Green's function given a pair of points
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual DerivativeTraits operator()(DerivativeTraits * source, DerivativeTraits * probe) const
    {
        // The distance has to be calculated using epsilonInv_ as metric:
        DerivativeTraits scratch = 0.0;
        Eigen::Matrix3d epsilonInv_ = this->profile_.epsilonInv();
        double detEps_ = this->profile_.detEps();
        for (int i = 0; i < 3; ++i) {
	        for (int j = 0; j < 3; ++j) {
		       scratch += (source[i] - probe[i]) * epsilonInv_(i, j) * (source[j] - probe[j]);
	        }
        }
        DerivativeTraits distance = sqrt(scratch);

        return (1.0/(sqrt(detEps_) * distance));
    }
    virtual std::ostream & printObject(std::ostream & os)
    {
        os << "Green's function type: anisotropic liquid" << std::endl;
        /*
        os << "Permittivity tensor diagonal (lab frame)   = " << epsilonLab_.transpose() << std::endl;
        os << "Permittivity tensor (molecule-fixed frame) =\n" << epsilon_ << std::endl;
        os << "Euler angles (molecule-to-lab frame)       = " << eulerAngles_.transpose();
        */
        return os;
    }
};

/*! \file AnisotropicLiquid.hpp
 *  \class AnisotropicLiquid
 *  \brief Green's functions for anisotropic liquid, described by a tensorial permittivity
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \tparam IntegratorPolicy policy for the calculation of diagonal elements
 *
 *  Explicit specialization
 */
template <typename IntegratorPolicy>
class AnisotropicLiquid<Numerical, IntegratorPolicy> : public GreensFunction<Numerical, IntegratorPolicy, Anisotropic,
                                     AnisotropicLiquid<Numerical, IntegratorPolicy> >
{
public:
    /*! \param[in] eigen_eps eigenvalues of the permittivity tensors
     *  \param[in] euler_ang Euler angles in degrees
     */
    AnisotropicLiquid(const Eigen::Vector3d & eigen_eps, const Eigen::Vector3d & euler_ang) :
        GreensFunction<Numerical, IntegratorPolicy, Anisotropic,
                  AnisotropicLiquid<Numerical, IntegratorPolicy> >() { this->profile_ = Anisotropic(eigen_eps, euler_ang); }
    virtual ~AnisotropicLiquid() {}
    /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e. the value of the
     *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    virtual double function(const Eigen::Vector3d & source, const Eigen::Vector3d & probe) const
    {
        Numerical sp[3], pp[3], res;
        sp[0] = source(0); sp[1] = source(1); sp[2] = source(2);
        pp[0] = probe(0);  pp[1] = probe(1);  pp[2] = probe(2);
        res = this->operator()(sp, pp);
        return res;
    }
    /*! Returns value of the kernel for the calculation of the \f$\mathcal{D}\f$ integral operator
     *  for the pair of points p1, p2:
     *  \f$ [\boldsymbol{\varepsilon}\nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)]\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  To obtain the kernel of the \f$\mathcal{D}^\dagger\f$ operator call this methods with \f$\mathbf{p}_1\f$
     *  and \f$\mathbf{p}_2\f$ exchanged and with \f$\mathbf{n}_{\mathbf{p}_2} = \mathbf{n}_{\mathbf{p}_1}\f$
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivative(const Eigen::Vector3d & /* direction */,
                              const Eigen::Vector3d & /* p1 */, const Eigen::Vector3d & /* p2 */) const
    {
        // Need the full gradient to get the kernel of D and D^\dagger
        /*
        Eigen::Vector3d scratch = this->profile_.epsilon() * (this->gradientProbe(p1, p2));
        return scratch.dot(direction);
        */
        return 0.0;
    }
    /*! Returns value of the directional derivative of the
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_1}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the source point.
     *  \param[in] normal_p1 the normal vector to p1
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivativeSource(const Eigen::Vector3d & normal_p1,
                            const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
    {
        using namespace std::placeholders;
        return threePointStencil(std::bind(&AnisotropicLiquid<Numerical, IntegratorPolicy>::function, this, _1, _2),
                                p1, p2, normal_p1, this->delta_);
    }
    /*! Returns value of the directional derivative of the
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point.
     *  \param[in] normal_p2 the normal vector to p2
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivativeProbe(const Eigen::Vector3d & normal_p2,
                                   const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
    {
        using namespace std::placeholders;
        return threePointStencil(std::bind(&AnisotropicLiquid<Numerical, IntegratorPolicy>::function, this, _1, _2),
                                p2, p1, normal_p2, this->delta_);
    }

    /*!
     *  Calculates the diagonal elements of the S operator: \f$ S_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalS(const Element & e) const
    {
        return this->diagonal_.computeS(*this, e);
    }
    /*!
     *  Calculates the diagonal elements of the D operator: \f$ D_{ii} \f$
     *  \param[in] e i-th finite element
     */
    virtual double diagonalD(const Element & e) const
    {
        return this->diagonal_.computeD(*this, e);
    }

    virtual double epsilon() const { return 0.0; }

    friend std::ostream & operator<<(std::ostream & os, AnisotropicLiquid & gf) {
        return gf.printObject(os);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    /*!
     *  Evaluates the Green's function given a pair of points
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual Numerical operator()(Numerical * source, Numerical * probe) const
    {
        // The distance has to be calculated using epsilonInv_ as metric:
        Numerical scratch = 0.0;
        Eigen::Matrix3d epsilonInv_ = this->profile_.epsilonInv();
        double detEps_ = this->profile_.detEps();
        for (int i = 0; i < 3; ++i) {
	        for (int j = 0; j < 3; ++j) {
		       scratch += (source[i] - probe[i]) * epsilonInv_(i, j) * (source[j] - probe[j]);
	        }
        }
        Numerical distance = sqrt(scratch);

        return (1.0/(sqrt(detEps_) * distance));
    }
    virtual std::ostream & printObject(std::ostream & os)
    {
        os << "Green's function type: anisotropic liquid" << std::endl;
        /*
        os << "Permittivity tensor diagonal (lab frame)   = " << epsilonLab_.transpose() << std::endl;
        os << "Permittivity tensor (molecule-fixed frame) =\n" << epsilon_ << std::endl;
        os << "Euler angles (molecule-to-lab frame)       = " << eulerAngles_.transpose();
        */
        return os;
    }
};

/*
namespace
{
    struct buildAnisotropicLiquid {
        template <typename DerivativeType>
        IGreensFunction * operator()(const greenData & _data) {
            DiagonalIntegrator * integrator =
		    DiagonalIntegratorFactory::TheDiagonalIntegratorFactory().createDiagonalIntegrator(_data.integratorType);
            return new AnisotropicLiquid<DerivativeType>(_data.epsilonTensor, _data.eulerAngles, integrator);
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
*/


#endif // ANISOTROPICLIQUID_HPP
