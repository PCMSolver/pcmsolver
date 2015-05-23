#ifndef COLLOCATIONINTEGRATOR_HPP
#define COLLOCATIONINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Dense>

class Element;

#include "DerivativeTypes.hpp"
#include "AnisotropicLiquid.hpp"
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

template <typename DerivativeTraits,
          typename ProfilePolicy>
struct CollocationIntegrator
{
    CollocationIntegrator() : factor_(1.07) {}
    ~CollocationIntegrator() {}

    double computeS(const Vacuum<DerivativeTraits, CollocationIntegrator> & /* gf */, const Element & e) const {
        double area = e.area();
	    return (factor_ * std::sqrt(4 * M_PI / area));
    }
    double computeD(const Vacuum<DerivativeTraits, CollocationIntegrator> & /* gf */, const Element & e) const {
        double area = e.area();
	    double radius = e.sphere().radius();
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
    }

    double computeS(const UniformDielectric<DerivativeTraits, CollocationIntegrator> & gf, const Element & e) const {
        double epsInv = 1.0 / gf.epsilon();
	    double area = e.area();
	    return (factor_ * std::sqrt(4 * M_PI / area) * epsInv);
    }
    double computeD(const UniformDielectric<DerivativeTraits, CollocationIntegrator> & gf, const Element & e) const {
        double epsInv = 1.0 / gf.epsilon();
	    double area = e.area();
	    double radius = e.sphere().radius();
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius) * epsInv);
    }

    double computeS(const IonicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const Element & /* e */) const {
        return 1.0;
    }
    double computeD(const IonicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const Element & /* e */) const {
        return 1.0;
    }

    double computeS(const AnisotropicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const Element & /* e */) const {
        return 1.0;
    }
    double computeD(const AnisotropicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const Element & /* e */) const {
        return 1.0;
    }

    /// Scaling factor for the collocation formulas
    double factor_;
};

#endif // COLLOCATIONINTEGRATOR_HPP
