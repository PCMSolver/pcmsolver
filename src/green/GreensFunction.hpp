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

#ifndef GREENSFUNCTION_HPP
#define GREENSFUNCTION_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "IGreensFunction.hpp"
#include "utils/Stencils.hpp"

/*! \file GreensFunction.hpp
 *  \class GreensFunction
 *  \brief Templated interface for Green's functions
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2014
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 *  \tparam IntegratorPolicy policy for the calculation of the matrix represenation of S and D
 *  \tparam ProfilePolicy    dielectric profile type
 */

template <typename DerivativeTraits,
          typename IntegratorPolicy,
          typename ProfilePolicy,
          typename Derived>
class GreensFunction: public IGreensFunction
{
public:
    GreensFunction() : delta_(1.0e-04), integrator_(IntegratorPolicy()) {}
    GreensFunction(double f) : delta_(1.0e-04), integrator_(IntegratorPolicy(f)) {}
    virtual ~GreensFunction() {}
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
        DerivativeTraits t1[3], t2[3];
        t1[0] = p1(0); t1[1] = p1(1); t1[2] = p1(2);
        t1[0][1] = normal_p1(0); t1[1][1] = normal_p1(1); t1[2][1] = normal_p1(2);
        t2[0] = p2(0); t2[1] = p2(1); t2[2] = p2(2);
        return this->operator()(t1, t2)[1];
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
                                   const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const __final
    {
        DerivativeTraits t1[3], t2[3];
        t1[0] = p1(0); t1[1] = p1(1); t1[2] = p1(2);
        t2[0] = p2(0); t2[1] = p2(1); t2[2] = p2(2);
        t2[0][1] = normal_p2(0); t2[1][1] = normal_p2(1); t2[2][1] = normal_p2(2);
        return this->operator()(t1, t2)[1];
    }
    /*! Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the source point.
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    Eigen::Vector3d gradientSource(const Eigen::Vector3d & p1,
                                           const Eigen::Vector3d & p2) const
    {
        return (Eigen::Vector3d() << derivativeSource(Eigen::Vector3d::UnitX(), p1, p2),
                derivativeSource(Eigen::Vector3d::UnitY(), p1, p2),
                derivativeSource(Eigen::Vector3d::UnitZ(), p1, p2)).finished();
    }
    /*! Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the probe point.
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    Eigen::Vector3d gradientProbe(const Eigen::Vector3d & p1,
                                          const Eigen::Vector3d & p2) const
    {
        return (Eigen::Vector3d() << derivativeProbe(Eigen::Vector3d::UnitX(), p1, p2),
                derivativeProbe(Eigen::Vector3d::UnitY(), p1, p2),
                derivativeProbe(Eigen::Vector3d::UnitZ(), p1, p2)).finished();
    }

    /*! Whether the Green's function describes a uniform environment */
    virtual bool uniform() const __final __override { return profiles::uniform(this->profile_); }
    /*! Returns a dielectric permittivity profile */
    virtual Permittivity permittivity() const __final __override { return this->profile_; }

    friend std::ostream & operator<<(std::ostream & os, GreensFunction & gf) {
        return gf.printObject(os);
    }
protected:
    /*! Evaluates the Green's function given a pair of points
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual DerivativeTraits operator()(DerivativeTraits * source, DerivativeTraits * probe) const = 0;
    /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e. the value of the
     *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     *  \note Relies on the implementation of operator() in the subclasses and that is all subclasses
     *  need to implement. Thus this method is marked __final.
     */
    virtual double kernelS_impl(const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const __final __override
    {
        DerivativeTraits sp[3], pp[3];
        sp[0] = p1(0); sp[1] = p1(1); sp[2] = p1(2);
        pp[0] = p2(0); pp[1] = p2(1); pp[2] = p2(2);
        return this->operator()(sp, pp)[0];
    }
    virtual std::ostream & printObject(std::ostream & os) __override
    {
        os << "Green's Function" << std::endl;
        return os;
    }
    double delta_;
    IntegratorPolicy integrator_;
    ProfilePolicy profile_;
};

template <typename IntegratorPolicy,
          typename ProfilePolicy,
          typename Derived>
class GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy, Derived>: public IGreensFunction
{
public:
    GreensFunction() : delta_(1.0e-04), integrator_(IntegratorPolicy()) {}
    GreensFunction(double f) : delta_(1.0e-04), integrator_(IntegratorPolicy(f)) {}
    virtual ~GreensFunction() {}
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
        return threePointStencil(pcm::bind(&GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy, Derived>::kernelS, this, pcm::_1, pcm::_2),
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
                                   const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const __final
    {
        return threePointStencil(pcm::bind(&GreensFunction<Numerical, IntegratorPolicy, ProfilePolicy, Derived>::kernelS, this, pcm::_1, pcm::_2),
                                p2, p1, normal_p2, this->delta_);
    }
    /*! Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the source point.
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    Eigen::Vector3d gradientSource(const Eigen::Vector3d & p1,
                                           const Eigen::Vector3d & p2) const
    {
        return (Eigen::Vector3d() << derivativeSource(Eigen::Vector3d::UnitX(), p1, p2),
                derivativeSource(Eigen::Vector3d::UnitY(), p1, p2),
                derivativeSource(Eigen::Vector3d::UnitZ(), p1, p2)).finished();
    }
    /*! Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the probe point.
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    Eigen::Vector3d gradientProbe(const Eigen::Vector3d & p1,
                                          const Eigen::Vector3d & p2) const
    {
        return (Eigen::Vector3d() << derivativeProbe(Eigen::Vector3d::UnitX(), p1, p2),
                derivativeProbe(Eigen::Vector3d::UnitY(), p1, p2),
                derivativeProbe(Eigen::Vector3d::UnitZ(), p1, p2)).finished();
    }

    /*! Whether the Green's function describes a uniform environment */
    virtual bool uniform() const __final __override { return profiles::uniform(this->profile_); }
    /*! Returns a dielectric permittivity profile */
    virtual Permittivity permittivity() const __final __override { return this->profile_; }

    friend std::ostream & operator<<(std::ostream & os, GreensFunction & gf) {
        return gf.printObject(os);
    }
protected:
    /*! Evaluates the Green's function given a pair of points
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual Numerical operator()(Numerical * source, Numerical * probe) const = 0;
    /*! Returns value of the kernel of the \f$\mathcal{S}\f$ integral operator, i.e. the value of the
     *  Greens's function for the pair of points p1, p2: \f$ G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     *  \note Relies on the implementation of operator() in the subclasses and that is all subclasses
     *  need to implement. Thus this method is marked __final.
     */
    virtual double kernelS_impl(const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const __final __override
    {
        return this->operator()(const_cast<Numerical *>(p1.data()), const_cast<Numerical *>(p2.data()));
    }
    virtual std::ostream & printObject(std::ostream & os) __override
    {
        os << "Green's Function" << std::endl;
        return os;
    }
    double delta_;
    IntegratorPolicy integrator_;
    ProfilePolicy profile_;
};

#endif // GREENSFUNCTION_HPP
