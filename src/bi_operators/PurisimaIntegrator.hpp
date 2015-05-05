#ifndef PURISIMAINTEGRATOR_HPP
#define PURISIMAINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Dense>

class Element;

#include "DerivativeTypes.hpp"
#include "DiagonalIntegrator.hpp"
#include "DiagonalIntegratorFactory.hpp"
#include "AnisotropicLiquid.hpp"
#include "IonicLiquid.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
#include "TanhSphericalDiffuse.hpp"

/*! \file PurisimaIntegrator.hpp
 *  \class PurisimaIntegrator
 *  \brief Implementation of diagonal elements of S and D using Purisima's formula
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  Calculates the diagonal elements of S as:
 *  \f[
 *  	S_{ii} = factor_ * \sqrt{\frac{4\pi}{a_i}}
 *  \f]
 *  while the diagonal elements of D are:
 *  \f[
 *  	D_{ii} = -(2\pi + \sum_{j\neq i} D_{ij}a_j)\frac{1}{a_i}
 *  \f]
 */

class PurisimaIntegrator : public DiagonalIntegrator
{
public:
    PurisimaIntegrator() : factor_(1.07) {}
    virtual ~PurisimaIntegrator() {}

    virtual double computeS(const Vacuum<double> * gf, const Element & e) const;
    virtual double computeS(const Vacuum<AD_directional> * gf, const Element & e) const;
    virtual double computeS(const Vacuum<AD_gradient> * gf, const Element & e) const;
    virtual double computeS(const Vacuum<AD_hessian> * gf, const Element & e) const;

    virtual double computeD(const Vacuum<double> * gf, const Element & e) const;
    virtual double computeD(const Vacuum<AD_directional> * gf, const Element & e) const;
    virtual double computeD(const Vacuum<AD_gradient> * gf, const Element & e) const;
    virtual double computeD(const Vacuum<AD_hessian> * gf, const Element & e) const;

    virtual double computeS(const UniformDielectric<double> * gf, const Element & e) const;
    virtual double computeS(const UniformDielectric<AD_directional> * gf, const Element & e) const;
    virtual double computeS(const UniformDielectric<AD_gradient> * gf, const Element & e) const;
    virtual double computeS(const UniformDielectric<AD_hessian> * gf, const Element & e) const;

    virtual double computeD(const UniformDielectric<double> * gf, const Element & e) const;
    virtual double computeD(const UniformDielectric<AD_directional> * gf, const Element & e) const;
    virtual double computeD(const UniformDielectric<AD_gradient> * gf, const Element & e) const;
    virtual double computeD(const UniformDielectric<AD_hessian> * gf, const Element & e) const;

    virtual double computeS(const IonicLiquid<double> * gf, const Element & e) const;
    virtual double computeS(const IonicLiquid<AD_directional> * gf, const Element & e) const;
    virtual double computeS(const IonicLiquid<AD_gradient> * gf, const Element & e) const;
    virtual double computeS(const IonicLiquid<AD_hessian> * gf, const Element & e) const;

    virtual double computeD(const IonicLiquid<double> * gf, const Element & e) const;
    virtual double computeD(const IonicLiquid<AD_directional> * gf, const Element & e) const;
    virtual double computeD(const IonicLiquid<AD_gradient> * gf, const Element & e) const;
    virtual double computeD(const IonicLiquid<AD_hessian> * gf, const Element & e) const;

    virtual double computeS(const AnisotropicLiquid<double> * gf, const Element & e) const;
    virtual double computeS(const AnisotropicLiquid<AD_directional> * gf, const Element & e) const;
    virtual double computeS(const AnisotropicLiquid<AD_gradient> * gf, const Element & e) const;
    virtual double computeS(const AnisotropicLiquid<AD_hessian> * gf, const Element & e) const;

    virtual double computeD(const AnisotropicLiquid<double> * gf, const Element & e) const;
    virtual double computeD(const AnisotropicLiquid<AD_directional> * gf, const Element & e) const;
    virtual double computeD(const AnisotropicLiquid<AD_gradient> * gf, const Element & e) const;
    virtual double computeD(const AnisotropicLiquid<AD_hessian> * gf, const Element & e) const;

    virtual double computeS(const TanhSphericalDiffuse * gf, const Element & e) const;

    virtual double computeD(const TanhSphericalDiffuse * gf, const Element & e) const;
private:
    /// Scaling factor for the collocation formula for S
    double factor_;
};

namespace
{
    DiagonalIntegrator * createPurisimaIntegrator()
    {
        return new PurisimaIntegrator();
    }
    const std::string PURISIMA("PURISIMA");
    const bool registeredPurisimaIntegrator = DiagonalIntegratorFactory::TheDiagonalIntegratorFactory().registerDiagonalIntegrator(
                                         PURISIMA, createPurisimaIntegrator);
}

#endif // PURISIMAINTEGRATOR_HPP
