#ifndef NUMERICALINTEGRATOR_HPP
#define NUMERICALINTEGRATOR_HPP

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

/*! \file NumericalIntegrator.hpp
 *  \class NumericalIntegrator
 *  \brief Implementation of diagonal elements of S and D using numerical integration
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  Calculates the diagonal elements of S and D by collocation, using numerical
 *  integration.
 */

class NumericalIntegrator : public DiagonalIntegrator
{
public:
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
};
namespace
{
    DiagonalIntegrator * createNumericalIntegrator()
    {
        return new NumericalIntegrator();
    }
    const std::string NUMERICAL("NUMERICAL");
    const bool registeredNumericalIntegrator = DiagonalIntegratorFactory::TheDiagonalIntegratorFactory().registerDiagonalIntegrator(
                                         NUMERICAL, createNumericalIntegrator);
}

#endif // NUMERICALINTEGRATOR_HPP
