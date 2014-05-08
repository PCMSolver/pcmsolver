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
    virtual double computeS(const Vacuum<double> * gf, double area) const = 0;
    virtual double computeS(const Vacuum<AD_directional> * gf, double area) const = 0;
    virtual double computeS(const Vacuum<AD_gradient> * gf, double area) const = 0;
    virtual double computeS(const Vacuum<AD_hessian> * gf, double area) const = 0;

    virtual double computeD(const Vacuum<double> * gf, double area, double radius) const = 0;
    virtual double computeD(const Vacuum<AD_directional> * gf, double area, double radius) const = 0;
    virtual double computeD(const Vacuum<AD_gradient> * gf, double area, double radius) const = 0;
    virtual double computeD(const Vacuum<AD_hessian> * gf, double area, double radius) const = 0;

    virtual double computeS(const UniformDielectric<double> * gf, double area) const = 0;
    virtual double computeS(const UniformDielectric<AD_directional> * gf, double area) const = 0;
    virtual double computeS(const UniformDielectric<AD_gradient> * gf, double area) const = 0;
    virtual double computeS(const UniformDielectric<AD_hessian> * gf, double area) const = 0;

    virtual double computeD(const UniformDielectric<double> * gf, double area, double radius) const = 0;
    virtual double computeD(const UniformDielectric<AD_directional> * gf, double area, double radius) const = 0;
    virtual double computeD(const UniformDielectric<AD_gradient> * gf, double area, double radius) const = 0;
    virtual double computeD(const UniformDielectric<AD_hessian> * gf, double area, double radius) const = 0;

    virtual double computeS(const IonicLiquid<double> * gf, double area) const = 0;
    virtual double computeS(const IonicLiquid<AD_directional> * gf, double area) const = 0;
    virtual double computeS(const IonicLiquid<AD_gradient> * gf, double area) const = 0;
    virtual double computeS(const IonicLiquid<AD_hessian> * gf, double area) const = 0;

    virtual double computeD(const IonicLiquid<double> * gf, double area, double radius) const = 0;
    virtual double computeD(const IonicLiquid<AD_directional> * gf, double area, double radius) const = 0;
    virtual double computeD(const IonicLiquid<AD_gradient> * gf, double area, double radius) const = 0;
    virtual double computeD(const IonicLiquid<AD_hessian> * gf, double area, double radius) const = 0;
};

#endif // DIAGONALINTEGRATOR_HPP
