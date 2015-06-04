/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef PURISIMAINTEGRATOR_HPP
#define PURISIMAINTEGRATOR_HPP

#include <cmath>
#include <functional>
#include <iosfwd>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

#include "IntegratorHelperFunctions.hpp"
#include "Element.hpp"
#include "AnisotropicLiquid.hpp"
#include "IonicLiquid.hpp"
#include "SphericalDiffuse.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"

/*! \file PurisimaIntegrator.hpp
 *  \struct PurisimaIntegrator
 *  \brief Implementation of the single and double layer operators matrix representation using
 *  one-point collocation and Purisima's strategy for the diagonal of D
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  Calculates the diagonal elements of S as:
 *  \f[
 *  	S_{ii} = factor_ * \sqrt{\frac{4\pi}{a_i}}
 *  \f]
 *  while the diagonal elements of D are:
 *  \f[
 *  	D_{ii} = -\left(2\pi + \sum_{j\neq i}D_{ij}a_j \right)\frac{1}{a_i}
 *  \f]
 *  The original reference is \cite Purisima1995
 */

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;

struct PurisimaIntegrator
{
    PurisimaIntegrator() : factor_(1.07) {}
    ~PurisimaIntegrator() {}

    /**@{ Single and double layer potentials for a Vacuum Green's function by collocation */
    /*! \tparam DerivativeTraits how the derivative of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const Vacuum<DerivativeTraits, PurisimaIntegrator> & gf, const std::vector<Element> & e) const {
        auto f = this->factor_;
        auto diagonal = [f] (const Element & el) -> double { return (f * std::sqrt(4 * M_PI / el.area())); };
        auto kernelS = std::bind(&Vacuum<DerivativeTraits, PurisimaIntegrator>::kernelS, gf, _1, _2);
        return integrator::singleLayer(e, diagonal, kernelS);
    }
    /*! \tparam DerivativeTraits how the derivative of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const Vacuum<DerivativeTraits, PurisimaIntegrator> & gf, const std::vector<Element> & e) const {
        auto f = this->factor_;
        auto diagonal = [f] (const Element & el) -> double { return (-f * std::sqrt(M_PI/ el.area()) * (1.0 / el.sphere().radius())); };
        auto kernelD = std::bind(&Vacuum<DerivativeTraits, PurisimaIntegrator>::kernelD, gf, _1, _2, _3);
        return integrator::doubleLayer(e, diagonal, kernelD);
    }
    /**@}*/

    /**@{ Single and double layer potentials for a UniformDielectric Green's function by collocation */
    /*! \tparam DerivativeTraits how the derivative of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const UniformDielectric<DerivativeTraits, PurisimaIntegrator> & gf, const std::vector<Element> & e) const {
        auto f = this->factor_;
        auto epsInv = 1.0 / gf.epsilon();
        auto diagonal = [f, epsInv] (const Element & el) -> double { return (f * std::sqrt(4 * M_PI / el.area()) * epsInv); };
        auto kernelS = std::bind(&UniformDielectric<DerivativeTraits, PurisimaIntegrator>::kernelS, gf, _1, _2);
        return integrator::singleLayer(e, diagonal, kernelS);
    }
    /*! \tparam DerivativeTraits how the derivative of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const UniformDielectric<DerivativeTraits, PurisimaIntegrator> & gf, const std::vector<Element> & e) const {
        auto f = this->factor_;
        auto diagonal = [f] (const Element & el) -> double { return (-f * std::sqrt(M_PI/ el.area()) * (1.0 / el.sphere().radius())); };
        auto kernelD = std::bind(&UniformDielectric<DerivativeTraits, PurisimaIntegrator>::kernelD, gf, _1, _2, _3);
        return integrator::doubleLayer(e, diagonal, kernelD);
    }
    /**@}*/

    /**@{ Single and double layer potentials for a IonicLiquid Green's function by collocation */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const IonicLiquid<DerivativeTraits, PurisimaIntegrator> & /* gf */, const std::vector<Element> & /* e */) const {
        throw std::runtime_error("Not implemented yet!");
    }
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const IonicLiquid<DerivativeTraits, PurisimaIntegrator> & /* gf */, const std::vector<Element> & /* e */) const {
        throw std::runtime_error("Not implemented yet!");
    }
    /**@}*/

    /**@{ Single and double layer potentials for an AnisotropicLiquid Green's function by collocation */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const AnisotropicLiquid<DerivativeTraits, PurisimaIntegrator> & /* gf */, const std::vector<Element> & /* e */) const {
        throw std::runtime_error("Not implemented yet!");
    }
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const AnisotropicLiquid<DerivativeTraits, PurisimaIntegrator> & /* gf */, const std::vector<Element> & /* e */) const {
        throw std::runtime_error("Not implemented yet!");
    }
    /**@}*/

    /**@{ Single and double layer potentials for a SphericalDiffuse Green's function by collocation */
    /*! \tparam ProfilePolicy the permittivity profile for the diffuse interface
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename ProfilePolicy>
    Eigen::MatrixXd singleLayer(const SphericalDiffuse<PurisimaIntegrator, ProfilePolicy> & gf, const std::vector<Element> & e) const {
        // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
        size_t mat_size = e.size();
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(mat_size, mat_size);
        for (size_t i = 0; i < mat_size; ++i) {
            // Fill diagonal
            // Diagonal of S inside the cavity
            double Sii_I = factor_ * std::sqrt(4 * M_PI / e[i].area());
            // "Diagonal" of Coulomb singularity separation coefficient
            double coulomb_coeff = gf.coefficientCoulomb(e[i].center(), e[i].center());
            // "Diagonal" of the image Green's function
            double image = gf.imagePotential(e[i].center(), e[i].center());
            S(i, i) = Sii_I / coulomb_coeff + image;
            Eigen::Vector3d source = e[i].center();
            for (size_t j = 0; j < mat_size; ++j) {
                // Fill off-diagonal
                Eigen::Vector3d probe = e[j].center();
                if (i != j) S(i, j) = gf.kernelS(source, probe);
            }
        }
        return S;
    }
    /*! \tparam ProfilePolicy the permittivity profile for the diffuse interface
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename ProfilePolicy>
    Eigen::MatrixXd doubleLayer(const SphericalDiffuse<PurisimaIntegrator, ProfilePolicy> & gf, const std::vector<Element> & e) const {
        // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
        size_t mat_size = e.size();
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(mat_size, mat_size);
        for (size_t i = 0; i < mat_size; ++i) {
            // Fill diagonal
            double area = e[i].area();
            double radius = e[i].sphere().radius();
            // Diagonal of S inside the cavity
            double Sii_I = factor_ * std::sqrt(4 * M_PI / area);
            // Diagonal of D inside the cavity
            double Dii_I = -factor_ * std::sqrt(M_PI/ area) * (1.0 / radius);
            // "Diagonal" of Coulomb singularity separation coefficient
            double coulomb_coeff = gf.coefficientCoulomb(e[i].center(), e[i].center());
            // "Diagonal" of the directional derivative of the Coulomb singularity separation coefficient
            double coeff_grad = gf.coefficientCoulombDerivative(e[i].normal(), e[i].center(), e[i].center()) / std::pow(coulomb_coeff, 2);
            // "Diagonal" of the directional derivative of the image Green's function
            double image_grad = gf.imagePotentialDerivative(e[i].normal(), e[i].center(), e[i].center());

            double eps_r2 = 0.0;
            std::tie(eps_r2, std::ignore) = gf.epsilon(e[i].center());

            D(i, i) = eps_r2 * (Dii_I / coulomb_coeff - Sii_I * coeff_grad + image_grad);
            Eigen::Vector3d source = e[i].center();
            for (size_t j = 0; j < mat_size; ++j) {
                // Fill off-diagonal
                Eigen::Vector3d probe = e[j].center();
                Eigen::Vector3d probeNormal = e[j].normal();
                probeNormal.normalize();
                if (i != j) D(i, j) = gf.kernelD(probeNormal, source, probe);
            }
        }
        return D;
    }
    /**@}*/

    /// Scaling factor for the collocation formulas
    double factor_;
};

#endif // PURISIMAINTEGRATOR_HPP
