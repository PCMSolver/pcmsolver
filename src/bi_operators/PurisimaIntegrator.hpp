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

#ifndef PURISIMAINTEGRATOR_HPP
#define PURISIMAINTEGRATOR_HPP

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

/*! \file PurisimaIntegrator.hpp
 *  \struct PurisimaIntegrator
 *  \brief Implementation of the single and double layer operators matrix representation using
 *  one-point collocation and Purisima's strategy for the diagonal of D
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  Calculates the diagonal elements of S as:
 *  \f[
 *  	S_{ii} = factor * \sqrt{\frac{4\pi}{a_i}}
 *  \f]
 *  while the diagonal elements of D are:
 *  \f[
 *  	D_{ii} = -\left(2\pi + \sum_{j\neq i}D_{ij}a_j \right)\frac{1}{a_i}
 *  \f]
 *  Note that the diagonal elements of S are obtained as in PurisimaIntegrator, only
 *  the calculation of the diagonal elements of D is different.
 *  The original reference is \cite Purisima1995
 */

struct PurisimaIntegrator
{
    PurisimaIntegrator() : factor(1.07) {}
    PurisimaIntegrator(double f) : factor(f) {}
    ~PurisimaIntegrator() {}

    /**@{ Single and double layer potentials for a Vacuum Green's function by collocation */
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const Vacuum<DerivativeTraits, PurisimaIntegrator> & gf, const std::vector<Element> & e) const {
        return integrator::singleLayer(e,
                pcm::bind(integrator::SI, this->factor, 1.0, pcm::_1),
                gf.exportKernelS());
    }
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const Vacuum<DerivativeTraits, PurisimaIntegrator> & gf, const std::vector<Element> & e) const {
        // Obtain off-diagonal first
        Eigen::MatrixXd D = offDiagonalD(e, gf.exportKernelD());
        // Fill diagonal based on Purisima's formula
        Eigen::VectorXd D_diag = diagonalD(e, D);
        D.diagonal() = D_diag;
        return D;
    }
    /**@}*/

    /**@{ Single and double layer potentials for a UniformDielectric Green's function by collocation */
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const UniformDielectric<DerivativeTraits, PurisimaIntegrator> & gf, const std::vector<Element> & e) const {
        return integrator::singleLayer(e,
                pcm::bind(integrator::SI, this->factor, gf.epsilon(), pcm::_1),
                gf.exportKernelS());
    }
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const UniformDielectric<DerivativeTraits, PurisimaIntegrator> & gf, const std::vector<Element> & e) const {
        // Obtain off-diagonal first
        Eigen::MatrixXd D = offDiagonalD(e, gf.exportKernelD());
        // Fill diagonal based on Purisima's formula
        Eigen::VectorXd D_diag = diagonalD(e, D);
        D.diagonal() = D_diag;
        return D;
    }
    /**@}*/

    /**@{ Single and double layer potentials for a IonicLiquid Green's function by collocation */
    template <typename DerivativeTraits>
      Eigen::MatrixXd singleLayer(const IonicLiquid<DerivativeTraits, PurisimaIntegrator> & /* gf */, const std::vector<Element> & e) const {
        PCMSOLVER_ERROR("PurisimaIntegrator::singleLayer not implemented yet for IonicLiquid", BOOST_CURRENT_FUNCTION);
        return Eigen::MatrixXd::Zero(e.size(), e.size());
      }
    template <typename DerivativeTraits>
      Eigen::MatrixXd doubleLayer(const IonicLiquid<DerivativeTraits, PurisimaIntegrator> & /* gf */, const std::vector<Element> & e) const {
        PCMSOLVER_ERROR("PurisimaIntegrator::doubleLayer not implemented yet for IonicLiquid", BOOST_CURRENT_FUNCTION);
        return Eigen::MatrixXd::Zero(e.size(), e.size());
      }
    /**@}*/

    /**@{ Single and double layer potentials for an AnisotropicLiquid Green's function by collocation */
    template <typename DerivativeTraits>
      Eigen::MatrixXd singleLayer(const AnisotropicLiquid<DerivativeTraits, PurisimaIntegrator> & /* gf */, const std::vector<Element> & e) const {
        PCMSOLVER_ERROR("PurisimaIntegrator::singleLayer not implemented yet for AnisotropicLiquid", BOOST_CURRENT_FUNCTION);
        return Eigen::MatrixXd::Zero(e.size(), e.size());

      }
    template <typename DerivativeTraits>
      Eigen::MatrixXd doubleLayer(const AnisotropicLiquid<DerivativeTraits, PurisimaIntegrator> & /* gf */, const std::vector<Element> & e) const {
        PCMSOLVER_ERROR("PurisimaIntegrator::doubleLayer not implemented yet for AnisotropicLiquid", BOOST_CURRENT_FUNCTION);
        return Eigen::MatrixXd::Zero(e.size(), e.size());
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
    Eigen::MatrixXd doubleLayer(const SphericalDiffuse<PurisimaIntegrator, ProfilePolicy> & gf, const std::vector<Element> & e) const {
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
    /*! Returns off-diagonal elements of the matrix representation of the double layer operator by collocation
     *  \param[in] elements list of finite elements
     *  \param[in] kernD    function for the evaluation of the off-diagonal of D
     */
    Eigen::MatrixXd offDiagonalD(const std::vector<Element> & elements,
            const integrator::KernelD & kernD) const
    {
        PCMSolverIndex mat_size = elements.size();
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(mat_size, mat_size);
        for (PCMSolverIndex i = 0; i < mat_size; ++i) {
            Eigen::Vector3d source = elements[i].center();
            for (PCMSolverIndex j = 0; j < mat_size; ++j) {
                // Fill off-diagonal
                Eigen::Vector3d probe = elements[j].center();
                Eigen::Vector3d probeNormal = elements[j].normal();
                probeNormal.normalize();
                if (i != j) D(i, j) = kernD(probeNormal, source, probe);
            }
        }
        return D;
    }
    /*! Returns diagonal elements of the matrix representation of the double layer operator by collocation
     *  \param[in] elements list of finite elements
     *  \param[in] D         the matrix representation of the double layer operator
     *  \f[
     *  	D_{ii} = -\left(2\pi + \sum_{j\neq i}D_{ij}a_j \right)\frac{1}{a_i}
     *  \f]
     */
    Eigen::VectorXd diagonalD(const std::vector<Element> & elements, const Eigen::MatrixXd & D) const
    {
        PCMSolverIndex mat_size = elements.size();
        Eigen::VectorXd D_diag = Eigen::VectorXd::Zero(mat_size);
        for (PCMSolverIndex i = 0; i < mat_size; ++i) {
            double D_ii = 0.0;
            for (PCMSolverIndex j = 0; j < mat_size; ++j) {
                if (j != i) D_ii += D(i, j) * elements[j].area();
            }
            D_diag(i) = - (2 * M_PI + D_ii) / (elements[i].area());
        }
        return D_diag;
    }
};

#endif // PURISIMAINTEGRATOR_HPP
