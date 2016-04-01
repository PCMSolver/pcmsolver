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

#ifndef NUMERICALINTEGRATOR_HPP
#define NUMERICALINTEGRATOR_HPP

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

/*! \file NumericalIntegrator.hpp
 *  \struct NumericalIntegrator
 *  \brief Implementation of the single and double layer operators matrix representation using one-point collocation
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  Calculates the diagonal elements of S and D by collocation, using numerical integration.
 */

struct NumericalIntegrator
{
    NumericalIntegrator() : factor(1.07) {}
    NumericalIntegrator(double f) : factor(f) {}
    ~NumericalIntegrator() {}
    /**@{ Single and double layer potentials for a Vacuum Green's function by collocation: numerical integration of diagonal */
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const Vacuum<DerivativeTraits, NumericalIntegrator> & gf, const std::vector<Element> & e) const {
        integrator::KernelS kernelS = gf.exportKernelS();
        integrator::Diagonal diagS = pcm::bind(&integrator::integrateS<32, 16>, kernelS, pcm::_1);
        return integrator::singleLayer(e, diagS, kernelS);
    }
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const Vacuum<DerivativeTraits, NumericalIntegrator> & gf, const std::vector<Element> & e) const {
        integrator::KernelD kernelD = gf.exportKernelD();
        integrator::Diagonal diagD = pcm::bind(&integrator::integrateD<32, 16>, kernelD, pcm::_1);
        return integrator::doubleLayer(e, diagD, kernelD);
    }
    /**@}*/

    /**@{ Single and double layer potentials for a UniformDielectric Green's function by collocation: numerical integration of diagonal */
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const UniformDielectric<DerivativeTraits, NumericalIntegrator> & gf, const std::vector<Element> & e) const {
        integrator::KernelS kernelS = gf.exportKernelS();
        integrator::Diagonal diagS = pcm::bind(&integrator::integrateS<32, 16>, kernelS, pcm::_1);
        return integrator::singleLayer(e, diagS, kernelS);
    }
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const UniformDielectric<DerivativeTraits, NumericalIntegrator> & gf, const std::vector<Element> & e) const {
        integrator::KernelD kernelD = gf.exportKernelD();
        integrator::Diagonal diagD = pcm::bind(&integrator::integrateD<32, 16>, kernelD, pcm::_1);
        return integrator::doubleLayer(e, diagD, kernelD);
    }
    /**@}*/

    /**@{ Single and double layer potentials for a IonicLiquid Green's function by collocation: numerical integration of diagonal */
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const IonicLiquid<DerivativeTraits, NumericalIntegrator> & gf, const std::vector<Element> & e) const {
        integrator::KernelS kernelS = gf.exportKernelS();
        integrator::Diagonal diagS = pcm::bind(&integrator::integrateS<32, 16>, kernelS, pcm::_1);
        return integrator::singleLayer(e, diagS, kernelS);
    }
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const IonicLiquid<DerivativeTraits, NumericalIntegrator> & gf, const std::vector<Element> & e) const {
        integrator::KernelD kernelD = gf.exportKernelD();
        integrator::Diagonal diagD = pcm::bind(&integrator::integrateD<32, 16>, kernelD, pcm::_1);
        return integrator::doubleLayer(e, diagD, kernelD);
    }
    /**@}*/

    /**@{ Single and double layer potentials for a AnisotropicLiquid Green's function by collocation: numerical integration of diagonal */
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd singleLayer(const AnisotropicLiquid<DerivativeTraits, NumericalIntegrator> & gf, const std::vector<Element> & e) const {
        integrator::KernelS kernelS = gf.exportKernelS();
        integrator::Diagonal diagS = pcm::bind(&integrator::integrateS<32, 16>, kernelS, pcm::_1);
        return integrator::singleLayer(e, diagS, kernelS);
    }
    /*! \tparam DerivativeTraits how the derivatives of the Greens's function are calculated
     *  \param[in] gf Green's function
     *  \param[in] e  list of finite elements
     */
    template <typename DerivativeTraits>
    Eigen::MatrixXd doubleLayer(const AnisotropicLiquid<DerivativeTraits, NumericalIntegrator> & gf, const std::vector<Element> & e) const {
        integrator::KernelD kernelD = gf.exportKernelD();
        integrator::Diagonal diagD = pcm::bind(&integrator::integrateD<32, 16>, kernelD, pcm::_1);
        return integrator::doubleLayer(e, diagD, kernelD);
    }
    /**@}*/

    /**@{ Single and double layer potentials for a SphericalDiffuse Green's function by collocation: numerical integration of diagonal */
    template <typename ProfilePolicy>
    Eigen::MatrixXd singleLayer(const SphericalDiffuse<NumericalIntegrator, ProfilePolicy> & /* gf */, const std::vector<Element> & e) const {
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
    Eigen::MatrixXd doubleLayer(const SphericalDiffuse<NumericalIntegrator, ProfilePolicy> & /* gf */, const std::vector<Element> & e) const {
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
//      pcm::tie(eps_r2, pcm::ignore) = gf.epsilon(e.center());

//      return eps_r2 * (Dii_I / coulomb_coeff - Sii_I * coeff_grad + image_grad);
        return Eigen::MatrixXd::Zero(e.size(), e.size());
    }
    /**@}*/
    /// Scaling factor for the collocation formulas (unused here)
    double factor;
};

#endif // NUMERICALINTEGRATOR_HPP
