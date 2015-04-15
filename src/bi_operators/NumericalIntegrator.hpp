#ifndef NUMERICALINTEGRATOR_HPP
#define NUMERICALINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Dense>

#include "DerivativeTypes.hpp"
#include "DiagonalIntegrator.hpp"
#include "DiagonalIntegratorFactory.hpp"
#include "IonicLiquid.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
#include "TanhSphericalDiffuse.hpp"

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
    virtual double computeS(const Vacuum<double> * gf, double area) const;
    virtual double computeS(const Vacuum<AD_directional> * gf, double area) const;
    virtual double computeS(const Vacuum<AD_gradient> * gf, double area) const;
    virtual double computeS(const Vacuum<AD_hessian> * gf, double area) const;

    virtual double computeD(const Vacuum<double> * gf, double area, double radius) const;
    virtual double computeD(const Vacuum<AD_directional> * gf, double area, double radius) const;
    virtual double computeD(const Vacuum<AD_gradient> * gf, double area, double radius) const;
    virtual double computeD(const Vacuum<AD_hessian> * gf, double area, double radius) const;

    virtual double computeS(const UniformDielectric<double> * gf, double area) const;
    virtual double computeS(const UniformDielectric<AD_directional> * gf, double area) const;
    virtual double computeS(const UniformDielectric<AD_gradient> * gf, double area) const;
    virtual double computeS(const UniformDielectric<AD_hessian> * gf, double area) const;

    virtual double computeD(const UniformDielectric<double> * gf, double area, double radius) const;
    virtual double computeD(const UniformDielectric<AD_directional> * gf, double area, double radius) const;
    virtual double computeD(const UniformDielectric<AD_gradient> * gf, double area, double radius) const;
    virtual double computeD(const UniformDielectric<AD_hessian> * gf, double area, double radius) const;

    virtual double computeS(const IonicLiquid<double> * gf, double area) const;
    virtual double computeS(const IonicLiquid<AD_directional> * gf, double area) const;
    virtual double computeS(const IonicLiquid<AD_gradient> * gf, double area) const;
    virtual double computeS(const IonicLiquid<AD_hessian> * gf, double area) const;

    virtual double computeD(const IonicLiquid<double> * gf, double area, double radius) const;
    virtual double computeD(const IonicLiquid<AD_directional> * gf, double area, double radius) const;
    virtual double computeD(const IonicLiquid<AD_gradient> * gf, double area, double radius) const;
    virtual double computeD(const IonicLiquid<AD_hessian> * gf, double area, double radius) const;

    virtual double computeS(const TanhSphericalDiffuse * gf, double area) const;

    virtual double computeD(const TanhSphericalDiffuse * gf, double area, double radius) const;
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
