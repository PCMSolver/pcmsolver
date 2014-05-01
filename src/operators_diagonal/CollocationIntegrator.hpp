#ifndef COLLOCATIONINTEGRATOR_HPP
#define COLLOCATIONINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "DiagonalIntegrator.hpp"
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
 *  	S_{ii} = 1.07\sqrt{\frac{4\pi}{a_i}}
 *  \f]
 *  while the diagonal elements of D are:
 *  \f[
 *  	D_{ii} = -1.07\sqrt{\frac{\pi}{a_i}} \frac{1}{R_I}
 *  \f]
 */

class CollocationIntegrator : public DiagonalIntegrator
{
public:
    virtual void computeS(const Vacuum<double> * gf) const = 0;
    virtual void computeS(const Vacuum<AD_directional> * gf) const = 0;
    virtual void computeS(const Vacuum<AD_gradient> * gf) const = 0;
    virtual void computeS(const Vacuum<AD_hessian> * gf) const = 0;

    virtual void computeD(const Vacuum<double> * gf) const = 0;
    virtual void computeD(const Vacuum<AD_directional> * gf) const = 0;
    virtual void computeD(const Vacuum<AD_gradient> * gf) const = 0;
    virtual void computeD(const Vacuum<AD_hessian> * gf) const = 0;

    virtual void computeS(const UniformDielectric<double> * gf) const = 0;
    virtual void computeS(const UniformDielectric<AD_directional> * gf) const = 0;
    virtual void computeS(const UniformDielectric<AD_gradient> * gf) const = 0;
    virtual void computeS(const UniformDielectric<AD_hessian> * gf) const = 0;

    virtual void computeD(const UniformDielectric<double> * gf) const = 0;
    virtual void computeD(const UniformDielectric<AD_directional> * gf) const = 0;
    virtual void computeD(const UniformDielectric<AD_gradient> * gf) const = 0;
    virtual void computeD(const UniformDielectric<AD_hessian> * gf) const = 0;

    virtual void computeS(const IonicLiquid<double> * gf) const = 0;
    virtual void computeS(const IonicLiquid<AD_directional> * gf) const = 0;
    virtual void computeS(const IonicLiquid<AD_gradient> * gf) const = 0;
    virtual void computeS(const IonicLiquid<AD_hessian> * gf) const = 0;

    virtual void computeD(const IonicLiquid<double> * gf) const = 0;
    virtual void computeD(const IonicLiquid<AD_directional> * gf) const = 0;
    virtual void computeD(const IonicLiquid<AD_gradient> * gf) const = 0;
    virtual void computeD(const IonicLiquid<AD_hessian> * gf) const = 0;
};

#endif // COLLOCATIONINTEGRATOR_HPP
