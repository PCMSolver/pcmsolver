/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#ifndef COLLOCATIONINTEGRATOR_HPP
#define COLLOCATIONINTEGRATOR_HPP

#include <cmath>
#include <iosfwd>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "IntegratorHelperFunctions.hpp"
#include "cavity/Element.hpp"
#include "green/AnisotropicLiquid.hpp"
#include "green/IonicLiquid.hpp"
#include "green/SphericalDiffuse.hpp"
#include "green/UniformDielectric.hpp"
#include "green/Vacuum.hpp"

/*! \file CollocationIntegrator.hpp
 *  \struct CollocationIntegrator
 *  \brief Implementation of the single and double layer operators matrix representation using one-point collocation
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  Calculates the diagonal elements of S as:
 *  \f[
 *  	S_{ii} = factor * \sqrt{\frac{4\pi}{a_i}}
 *  \f]
 *  while the diagonal elements of D are:
 *  \f[
 *  	D_{ii} = -factor * \sqrt{\frac{\pi}{a_i}} \frac{1}{R_I}
 *  \f]
 */

struct CollocationIntegrator
{
    CollocationIntegrator() : factor(1.07) {}
    CollocationIntegrator(double f) : factor(f) {}
    ~CollocationIntegrator() {}

    /**@{ Single and double layer potentials for a Vacuum Green's function by collocation */
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const Vacuum<DerivativeTraits, CollocationIntegrator> & gf, const std::vector<Element> & e) const {
        return integrator::singleLayer(e,
                pcm::bind(integrator::SI, this->factor, 1.0, pcm::_1),
                gf.exportKernelS());
    }
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const Vacuum<DerivativeTraits, CollocationIntegrator> & gf, const std::vector<Element> & e) const {
        return integrator::doubleLayer(e,
                                       pcm::bind(integrator::DI, this->factor, pcm::_1),
                                       gf.exportKernelD());
    }
    /**@}*/

    /**@{ Single and double layer potentials for a UniformDielectric Green's function by collocation */
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const UniformDielectric<DerivativeTraits, CollocationIntegrator> & gf, const std::vector<Element> & e) const {
        return integrator::singleLayer(e,
                pcm::bind(integrator::SI, this->factor, gf.epsilon(), pcm::_1),
                gf.exportKernelS());
    }
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const UniformDielectric<DerivativeTraits, CollocationIntegrator> & gf, const std::vector<Element> & e) const {
        return integrator::doubleLayer(e,
                                       pcm::bind(integrator::DI, this->factor, pcm::_1),
                                       gf.exportKernelD());
    }
    /**@}*/

    /**@{ Single and double layer potentials for a IonicLiquid Green's function by collocation */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const IonicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const std::vector<Element> & e) const {
        PCMSOLVER_ERROR("CollocationIntegrator::singleLayer not implemented yet for IonicLiquid", BOOST_CURRENT_FUNCTION);
        return Eigen::MatrixXd::Zero(e.size(), e.size());
    }
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const IonicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const std::vector<Element> & e) const {
        PCMSOLVER_ERROR("CollocationIntegrator::doubleLayer not implemented yet for IonicLiquid", BOOST_CURRENT_FUNCTION);
        return Eigen::MatrixXd::Zero(e.size(), e.size());
    }
    /**@}*/

    /**@{ Single and double layer potentials for an AnisotropicLiquid Green's function by collocation */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const AnisotropicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const std::vector<Element> & e) const {
        PCMSOLVER_ERROR("CollocationIntegrator::singleLayer not implemented yet for AnisotropicLiquid", BOOST_CURRENT_FUNCTION);
        return Eigen::MatrixXd::Zero(e.size(), e.size());
    }
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const AnisotropicLiquid<DerivativeTraits, CollocationIntegrator> & /* gf */, const std::vector<Element> & e) const {
        PCMSOLVER_ERROR("CollocationIntegrator::doubleLayer not implemented yet for AnisotropicLiquid", BOOST_CURRENT_FUNCTION);
        return Eigen::MatrixXd::Zero(e.size(), e.size());
    }
    /**@}*/

    /**@{ Single and double layer potentials for a SphericalDiffuse Green's function by collocation */
    /*! \tparam ProfilePolicy the permittivity profile for the diffuse interface
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename ProfilePolicy>
    Eigen::MatrixXd singleLayer(const SphericalDiffuse<CollocationIntegrator, ProfilePolicy> & gf, const std::vector<Element> & e) const {
        // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
        PCMSolverIndex mat_size = e.size();
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(mat_size, mat_size);
        for (PCMSolverIndex i = 0; i < mat_size; ++i) {
            // Fill diagonal
            // Diagonal of S inside the cavity
            double Sii_I = factor * std::sqrt(4 * M_PI / e[i].area());
            // "Diagonal" of Coulomb singularity separation coefficient
            double coulomb_coeff = gf.coefficientCoulomb(e[i].center(), e[i].center());
            // "Diagonal" of the image Green's function
            double image = gf.imagePotential(e[i].center(), e[i].center());
            S(i, i) = Sii_I / coulomb_coeff + image;
            Eigen::Vector3d source = e[i].center();
            for (PCMSolverIndex j = 0; j < mat_size; ++j) {
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
    Eigen::MatrixXd doubleLayer(const SphericalDiffuse<CollocationIntegrator, ProfilePolicy> & gf, const std::vector<Element> & e) const {
        // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
        PCMSolverIndex mat_size = e.size();
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(mat_size, mat_size);
        for (PCMSolverIndex i = 0; i < mat_size; ++i) {
            // Fill diagonal
            double area = e[i].area();
            double radius = e[i].sphere().radius;
            // Diagonal of S inside the cavity
            double Sii_I = factor * std::sqrt(4 * M_PI / area);
            // Diagonal of D inside the cavity
            double Dii_I = -factor * std::sqrt(M_PI/ area) * (1.0 / radius);
            // "Diagonal" of Coulomb singularity separation coefficient
            double coulomb_coeff = gf.coefficientCoulomb(e[i].center(), e[i].center());
            // "Diagonal" of the directional derivative of the Coulomb singularity separation coefficient
            double coeff_grad = gf.coefficientCoulombDerivative(e[i].normal(), e[i].center(), e[i].center()) / std::pow(coulomb_coeff, 2);
            // "Diagonal" of the directional derivative of the image Green's function
            double image_grad = gf.imagePotentialDerivative(e[i].normal(), e[i].center(), e[i].center());

            double eps_r2 = 0.0;
            pcm::tie(eps_r2, pcm::ignore) = gf.epsilon(e[i].center());

            D(i, i) = eps_r2 * (Dii_I / coulomb_coeff - Sii_I * coeff_grad + image_grad);
            Eigen::Vector3d source = e[i].center();
            for (PCMSolverIndex j = 0; j < mat_size; ++j) {
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
    double factor;
};

#endif // COLLOCATIONINTEGRATOR_HPP
