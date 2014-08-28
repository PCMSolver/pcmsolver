#ifndef COLLOCATIONINTEGRATOR_HPP
#define COLLOCATIONINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Dense>

#include "DerivativeTypes.hpp"
#include "DiagonalIntegrator.hpp"
#include "DiagonalIntegratorFactory.hpp"
#include "IonicLiquid.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"

/*! \file CollocationIntegrator.hpp
 *  \class CollocationIntegrator
 *  \brief Implementation of diagonal elements of S and D using approximate collocation
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  Calculates the diagonal elements of S as:
 *  \f[
 *  	S_{ii} = factor_ * \sqrt{\frac{4\pi}{a_i}}
 *  \f]
 *  while the diagonal elements of D are:
 *  \f[
 *  	D_{ii} = -factor_ * \sqrt{\frac{\pi}{a_i}} \frac{1}{R_I}
 *  \f]
 */

class CollocationIntegrator : public DiagonalIntegrator
{
public:
    CollocationIntegrator() : factor_(1.07) {}
    virtual ~CollocationIntegrator() {}

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
private:
    double factor_;
};

namespace
{
    DiagonalIntegrator * createCollocationIntegrator()
    {
        return new CollocationIntegrator();
    }
    const std::string COLLOCATION("COLLOCATION");
    const bool registeredCollocationIntegrator = DiagonalIntegratorFactory::TheDiagonalIntegratorFactory().registerDiagonalIntegrator(
                                         COLLOCATION, createCollocationIntegrator);
}

#endif // COLLOCATIONINTEGRATOR_HPP
