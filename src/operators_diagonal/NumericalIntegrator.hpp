#ifndef NUMERICALINTEGRATOR_HPP
#define NUMERICALINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "DiagonalIntegrator.hpp"
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

#endif // NUMERICALINTEGRATOR_HPP
