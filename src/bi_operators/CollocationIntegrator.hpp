#ifndef COLLOCATIONINTEGRATOR_HPP
#define COLLOCATIONINTEGRATOR_HPP

#include <cmath>
#include <functional>
#include <iosfwd>
#include <stdexcept>

#include "Config.hpp"

#include <Eigen/Dense>

#include "IntegratorHelperFunctions.hpp"
#include "Element.hpp"
#include "AnisotropicLiquid.hpp"
#include "IonicLiquid.hpp"
#include "SphericalDiffuse.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"

/*! \file CollocationIntegrator.hpp
 *  \struct CollocationIntegrator
 *  \brief Implementation of the single and double layer operators matrix representation using one-point collocation
 *  \author Roberto Di Remigio
 *  \date 2015
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

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;
using integrator::KernelS;
using integrator::KernelD;

struct CollocationIntegrator
{
    CollocationIntegrator() : factor_(1.07) {}
    ~CollocationIntegrator() {}

    /**@{ Single and double layer potentials for a Vacuum Green's function by collocation */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const Vacuum<DerivativeTraits, CollocationIntegrator> & gf, const std::vector<Element> & e) const {
        auto f = this->factor_;
        auto diagonal = [f] (double area) -> double { return (f * std::sqrt(4 * M_PI / area)); };
        KernelS kernelS = std::bind(&Vacuum<DerivativeTraits, CollocationIntegrator>::kernelS, gf, _1, _2);
        return integrator::singleLayer(e, diagonal, kernelS);
    }
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const Vacuum<DerivativeTraits, CollocationIntegrator> & gf, const std::vector<Element> & e) const {
        auto f = this->factor_;
        auto diagonal = [f] (double area, double radius) -> double { return (-f * std::sqrt(M_PI/ area) * (1.0 / radius)); };
        KernelD kernelD = std::bind(&Vacuum<DerivativeTraits, CollocationIntegrator>::kernelD, gf, _1, _2, _3);
        return integrator::doubleLayer(e, diagonal, kernelD);
    }
    /**@}*/

    /**@{ Single and double layer potentials for a UniformDielectric Green's function by collocation */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const UniformDielectric<DerivativeTraits, CollocationIntegrator> & gf, const std::vector<Element> & e) const {
        auto f = this->factor_;
        auto epsInv = 1.0 / gf.epsilon();
        auto diagonal = [f, epsInv] (double area) -> double { return (f * std::sqrt(4 * M_PI / area) * epsInv); };
        KernelS kernelS = std::bind(&UniformDielectric<DerivativeTraits, CollocationIntegrator>::kernelS, gf, _1, _2);
        return integrator::singleLayer(e, diagonal, kernelS);
    }
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const UniformDielectric<DerivativeTraits, CollocationIntegrator> & gf, const std::vector<Element> & e) const {
        auto f = this->factor_;
        auto diagonal = [f] (double area, double radius) -> double { return (-f * std::sqrt(M_PI/ area) * (1.0 / radius)); };
        KernelD kernelD = std::bind(&UniformDielectric<DerivativeTraits, CollocationIntegrator>::kernelD, gf, _1, _2, _3);
        return integrator::doubleLayer(e, diagonal, kernelD);
    }
    /**@}*/

    /**@{ Single and double layer potentials for a IonicLiquid Green's function by collocation */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const IonicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const std::vector<Element> & /* e */) const {
        throw std::runtime_error("Not implemented yet!");
    }
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const IonicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const std::vector<Element> & /* e */) const {
        throw std::runtime_error("Not implemented yet!");
    }
    /**@}*/

    /**@{ Single and double layer potentials for an AnisotropicLiquid Green's function by collocation */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const AnisotropicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const std::vector<Element> & /* e */) const {
        throw std::runtime_error("Not implemented yet!");
    }
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const AnisotropicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const std::vector<Element> & /* e */) const {
        throw std::runtime_error("Not implemented yet!");
    }
    /**@}*/

    /**@{ Single and double layer potentials for a SphericalDiffuse Green's function by collocation */
    template <typename ProfilePolicy>
    Eigen::MatrixXd singleLayer(const SphericalDiffuse<CollocationIntegrator, ProfilePolicy> & /* gf */, const std::vector<Element> & e) const {
//      // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
//      double area = e.area();
//      // Diagonal of S inside the cavity
//      double Sii_I = factor_ * std::sqrt(4 * M_PI / area);
//      // "Diagonal" of Coulomb singularity separation coefficient
//      double coulomb_coeff = gf.coefficientCoulomb(e.center(), e.center());
//      // "Diagonal" of the image Green's function
//      double image = gf.imagePotential(e.center(), e.center());

//      return (Sii_I / coulomb_coeff + image);
        return Eigen::MatrixXd::Zero(e.size(), e.size());
    }
    template <typename ProfilePolicy>
    Eigen::MatrixXd doubleLayer(const SphericalDiffuse<CollocationIntegrator, ProfilePolicy> & /* gf */, const std::vector<Element> & e) const {
//      // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
//      double area = e.area();
//      double radius = e.sphere().radius();
//      // Diagonal of S inside the cavity
//      double Sii_I = factor_ * std::sqrt(4 * M_PI / area);
//      // Diagonal of D inside the cavity
//      double Dii_I = -factor_ * std::sqrt(M_PI/ area) * (1.0 / radius);
//      // "Diagonal" of Coulomb singularity separation coefficient
//      double coulomb_coeff = gf.coefficientCoulomb(e.center(), e.center());
//      // "Diagonal" of the directional derivative of the Coulomb singularity separation coefficient
//      double coeff_grad = gf.coefficientCoulombDerivative(e.normal(), e.center(), e.center()) / std::pow(coulomb_coeff, 2);
//      // "Diagonal" of the directional derivative of the image Green's function
//      double image_grad = gf.imagePotentialDerivative(e.normal(), e.center(), e.center());

//      double eps_r2 = 0.0;
//      std::tie(eps_r2, std::ignore) = gf.epsilon(e.center());

//      return eps_r2 * (Dii_I / coulomb_coeff - Sii_I * coeff_grad + image_grad);
        return Eigen::MatrixXd::Zero(e.size(), e.size());
    }
    /**@}*/

    /// Scaling factor for the collocation formulas
    double factor_;
};

#endif // COLLOCATIONINTEGRATOR_HPP
