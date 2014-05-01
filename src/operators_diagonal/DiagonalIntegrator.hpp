#ifndef DIAGONALINTEGRATOR_HPP
#define DIAGONALINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "IonicLiquid.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"

/*! \file DiagonalIntegrator.hpp
 *  \class DiagonalIntegrator
 *  \brief Abstract Base Class for implementation of diagonal elements of S and D
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  This class encapsulates the calculation of the diagonal elements of the S and D
 *  operators needed to set up the PCM equations.
 *  Based on the ideas of the Strategy Pattern.
 */

class DiagonalIntegrator
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

#endif // DIAGONALINTEGRATOR_HPP
