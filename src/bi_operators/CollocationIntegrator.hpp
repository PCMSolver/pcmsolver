#ifndef COLLOCATIONINTEGRATOR_HPP
#define COLLOCATIONINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Dense>

class Element;

#include "DerivativeTypes.hpp"
#include "AnisotropicLiquid.hpp"
#include "IonicLiquid.hpp"
#include "SphericalDiffuse.hpp"
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

    double computeS(const SphericalDiffuse<CollocationIntegrator, ProfilePolicy> & gf, const Element & e) const {
        // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
        double area = e.area();
        // Diagonal of S inside the cavity
        double Sii_I = factor_ * std::sqrt(4 * M_PI / area);
        // "Diagonal" of Coulomb singularity separation coefficient
        double coulomb_coeff = gf.coefficientCoulomb(e.center(), e.center());
        // "Diagonal" of the image Green's function
        double image = gf.imagePotential(e.center(), e.center());

        return (Sii_I / coulomb_coeff + image);
    }
    double computeD(const SphericalDiffuse<CollocationIntegrator, ProfilePolicy> & gf, const Element & e) const {
        // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
        double area = e.area();
        double radius = e.sphere().radius();
        // Diagonal of S inside the cavity
        double Sii_I = factor_ * std::sqrt(4 * M_PI / area);
        // Diagonal of D inside the cavity
        double Dii_I = -factor_ * std::sqrt(M_PI/ area) * (1.0 / radius);
        // "Diagonal" of Coulomb singularity separation coefficient
        double coulomb_coeff = gf.coefficientCoulomb(e.center(), e.center());
        // "Diagonal" of the directional derivative of the Coulomb singularity separation coefficient
        double coeff_grad = gf.coefficientCoulombDerivative(e.normal(), e.center(), e.center()) / std::pow(coulomb_coeff, 2);
        // "Diagonal" of the directional derivative of the image Green's function
        double image_grad = gf.imagePotentialDerivative(e.normal(), e.center(), e.center());

        return (Dii_I / coulomb_coeff - Sii_I * coeff_grad + image_grad);
    }

    /// Scaling factor for the collocation formulas
    double factor_;
};

#endif // COLLOCATIONINTEGRATOR_HPP
