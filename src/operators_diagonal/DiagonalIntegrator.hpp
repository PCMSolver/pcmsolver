#ifndef DIAGONALINTEGRATOR_HPP
#define DIAGONALINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "Element.hpp"
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
    virtual double computeS(const Vacuum<double> * gf, const Element & e) const = 0;
    virtual double computeS(const Vacuum<AD_directional> * gf, const Element & e) const = 0;
    virtual double computeS(const Vacuum<AD_gradient> * gf, const Element & e) const = 0;
    virtual double computeS(const Vacuum<AD_hessian> * gf, const Element & e) const = 0;

    virtual double computeD(const Vacuum<double> * gf, const Element & e) const = 0;
    virtual double computeD(const Vacuum<AD_directional> * gf, const Element & e) const = 0;
    virtual double computeD(const Vacuum<AD_gradient> * gf, const Element & e) const = 0;
    virtual double computeD(const Vacuum<AD_hessian> * gf, const Element & e) const = 0;

    virtual double computeS(const UniformDielectric<double> * gf, const Element & e) const = 0;
    virtual double computeS(const UniformDielectric<AD_directional> * gf, const Element & e) const = 0;
    virtual double computeS(const UniformDielectric<AD_gradient> * gf, const Element & e) const = 0;
    virtual double computeS(const UniformDielectric<AD_hessian> * gf, const Element & e) const = 0;

    virtual double computeD(const UniformDielectric<double> * gf, const Element & e) const = 0;
    virtual double computeD(const UniformDielectric<AD_directional> * gf, const Element & e) const = 0;
    virtual double computeD(const UniformDielectric<AD_gradient> * gf, const Element & e) const = 0;
    virtual double computeD(const UniformDielectric<AD_hessian> * gf, const Element & e) const = 0;

    virtual double computeS(const IonicLiquid<double> * gf, const Element & e) const = 0;
    virtual double computeS(const IonicLiquid<AD_directional> * gf, const Element & e) const = 0;
    virtual double computeS(const IonicLiquid<AD_gradient> * gf, const Element & e) const = 0;
    virtual double computeS(const IonicLiquid<AD_hessian> * gf, const Element & e) const = 0;

    virtual double computeD(const IonicLiquid<double> * gf, const Element & e) const = 0;
    virtual double computeD(const IonicLiquid<AD_directional> * gf, const Element & e) const = 0;
    virtual double computeD(const IonicLiquid<AD_gradient> * gf, const Element & e) const = 0;
    virtual double computeD(const IonicLiquid<AD_hessian> * gf, const Element & e) const = 0;
};

#endif // DIAGONALINTEGRATOR_HPP
