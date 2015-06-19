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
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#ifndef ANALYTICEVALUATE_HPP
#define ANALYTICEVALUATE_HPP

#include <cmath>

#include "Config.hpp"

#include <Eigen/Core>

#include "MathUtils.hpp"

/*! \brief Analytic evaluation of vacuum Green's function and its derivatives
 *  
 *  \f[
 *  \begin{align}
 *     G(\vect{r},\vect{r}^\prime) &= \frac{1}{|\vect{r}-\vect{r}^\prime|} \\
 *     \pderiv{}{{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) &= \frac{(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}^\prime}}{|\vect{r}-\vect{r}^\prime|^3} \\
 *     \pderiv{}{{\vect{n}_{\vect{r}}}}G(\vect{r},\vect{r}^\prime) &= -\frac{(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}}}{|\vect{r}-\vect{r}^\prime|^3} \\
 *     \frac{\partial^2}{\partial{\vect{n}_{\vect{r}}}\partial{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) &=
 *     \frac{\vect{n}_{\vect{r}}\cdot \vect{n}_{\vect{r}^\prime}}{|\vect{r}-\vect{r}^\prime|^3}
 *     -3\frac{[(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}^\prime}]}{|\vect{r}-\vect{r}^\prime|^5}
 *  \end{align}
 *  \f]	      
 */
inline Eigen::Array4d analyticVacuum(const Eigen::Vector3d & spNormal,
                                const Eigen::Vector3d & sp,
                                const Eigen::Vector3d & ppNormal, const Eigen::Vector3d & pp) {
    Eigen::Array4d result = Eigen::Array4d::Zero();
    double distance = (sp - pp).norm();
    double distance_3 = std::pow(distance, 3);
    double distance_5 = std::pow(distance, 5);

    // Value of the function
    result(0) = 1.0 / distance;
    // Value of the directional derivative wrt probe
    result(1) = (sp - pp).dot(ppNormal) / distance_3 ;
    // Directional derivative wrt source
    result(2) = - (sp - pp).dot(spNormal) / distance_3;
    // Value of the Hessian
    result(3) = spNormal.dot(ppNormal) / distance_3 - 3 * ((
                    sp - pp).dot(spNormal))*((sp - pp).dot(
                                ppNormal)) / distance_5;

    return result;
}

/*! \brief Analytic evaluation of uniform dielectric Green's function and its derivatives
 *
 *  \f[
 *  \begin{align}
 *     G(\vect{r},\vect{r}^\prime) &= \frac{1}{|\diel\vect{r}-\vect{r}^\prime|} \\
 *     \pderiv{}{{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) &= \frac{(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}^\prime}}{\diel|\vect{r}-\vect{r}^\prime|^3} \\
 *     \pderiv{}{{\vect{n}_{\vect{r}}}}G(\vect{r},\vect{r}^\prime) &= -\frac{(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}}}{\diel|\vect{r}-\vect{r}^\prime|^3} \\
 *     \frac{\partial^2}{\partial{\vect{n}_{\vect{r}}}\partial{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) &=
 *     \frac{\vect{n}_{\vect{r}}\cdot \vect{n}_{\vect{r}^\prime}}{\diel|\vect{r}-\vect{r}^\prime|^3}
 *     -3\frac{[(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}^\prime}]}{\diel|\vect{r}-\vect{r}^\prime|^5}
 *  \end{align}
 *  \f]	      
 */
inline Eigen::Array4d analyticUniformDielectric(double eps, const Eigen::Vector3d & spNormal,
                                const Eigen::Vector3d & sp,
                                const Eigen::Vector3d & ppNormal, const Eigen::Vector3d & pp) {
    Eigen::Array4d result = Eigen::Array4d::Zero();
    double distance = (sp - pp).norm();
    double distance_3 = std::pow(distance, 3);
    double distance_5 = std::pow(distance, 5);

    // Value of the function
    result(0) = 1.0 / (eps * distance);
    // Value of the directional derivative wrt probe
    result(1) = (sp - pp).dot(ppNormal) / (eps * distance_3);
    // Directional derivative wrt source
    result(2) = - (sp - pp).dot(spNormal) / (eps * distance_3);
    // Value of the Hessian
    result(3) = spNormal.dot(ppNormal) / (eps * distance_3) - 3 * ((
                    sp - pp).dot(spNormal))*((sp - pp).dot(
                                ppNormal)) / (eps * distance_5);

    return result;
}

/*! \brief Analytic evaluation of ionic liquid Green's function and its derivatives
 *
 *  \f[
 *  \begin{align}
 *  \phantom{G(\vect{r},\vect{r}^\prime)}
 *  &\begin{aligned}
 *    G(\vect{r},\vect{r}^\prime) = \frac{\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|}
 *  \end{aligned}\\
 *  &\begin{aligned}  
 *    \pderiv{}{{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) = 
 *    \frac{(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}^\prime}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^3} 
 *    +\kappa\frac{(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}^\prime}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^2}
 *   \end{aligned}\\
 *   &\begin{aligned}  
 *     \pderiv{}{{\vect{n}_{\vect{r}}}}G(\vect{r},\vect{r}^\prime) = 
 *     -\frac{(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^3} 
 *     -\kappa\frac{(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^2}
 *    \end{aligned}\\
 *   &\begin{aligned}
 *     \frac{\partial^2}{\partial{\vect{n}_{\vect{r}}}\partial{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) &=
 *     \frac{\vect{n}_{\vect{r}}\cdot \vect{n}_{\vect{r}^\prime}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^3}
 *     -\kappa\frac{[(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}^\prime}]\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^4}\\
 *     &-3\frac{[(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}^\prime}]\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^5}
 *     + \kappa\frac{\vect{n}_{\vect{r}}\cdot \vect{n}_{\vect{r}^\prime}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^2} \\
 *     &-\kappa^2\frac{[(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}^\prime}]\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^3}\\
 *     &-2\kappa\frac{[(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot \vect{n}_{\vect{r}^\prime}]\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^4}
 *    \end{aligned}
 *  \end{align}
 *  \f]
 */
inline Eigen::Array4d analyticIonicLiquid(double eps, double k,
                                const Eigen::Vector3d & spNormal,
                                const Eigen::Vector3d & sp,
                                const Eigen::Vector3d & ppNormal, const Eigen::Vector3d & pp) {
    Eigen::Array4d result = Eigen::Array4d::Zero();
    double distance = (sp - pp).norm();
    double distance_3 = std::pow(distance, 3);
    double distance_5 = std::pow(distance, 5);

    // Value of the function
    result(0) = std::exp(- k * distance) / (eps * distance);
    // Value of the directional derivative wrt probe
    result(1) = (sp - pp).dot(ppNormal) * (1 + k * distance ) * std::exp(
                    - k * distance) / (eps * distance_3);
    // Directional derivative wrt source
    result(2) = - (sp - pp).dot(spNormal) * (1 + k * distance ) * std::exp(
                    - k * distance) / (eps * distance_3);
    // Value of the Hessian
    result(3) = spNormal.dot(ppNormal) * (1 + k * distance) * std::exp(
                    - k * distance) / (eps * distance_3)
                - std::pow(k, 2) * (sp - pp).dot(spNormal) * (sp - pp).dot(
                    ppNormal) * std::exp(- k * distance) / (eps * distance_3)
                - 3 * (sp - pp).dot(spNormal) * (sp - pp).dot(
                    ppNormal) * (1 + k * distance) * std::exp(- k * distance) /
                (eps * distance_5);

    return result;
}

/*! \brief Analytic evaluation of anisotropic liquid Green's function and its derivatives
 *
 *  \f[
 *  \begin{align}
 *   \phantom{G(\vect{r},\vect{r}^\prime)}
 *   &\begin{aligned}
 *     G(\vect{r},\vect{r}^\prime) = \frac{1}{4\pi\sqrt{\det{\bm{\diel}}}\sqrt{(\vect{r}-\vect{r}^\prime)^t\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)}} 
 *     \end{aligned}\\
 *   &\begin{aligned}  
 *     \pderiv{}{{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) = \frac{[\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]\cdot \vect{n}_{\vect{r}^\prime}}{4\pi\sqrt{\det{\bm{\diel}}}
 *     [(\vect{r}-\vect{r}^\prime)^t\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]^{\frac{3}{2}}} 
 *     \end{aligned}\\
 *   &\begin{aligned}
 *     \pderiv{}{{\vect{n}_{\vect{r}}}}G(\vect{r},\vect{r}^\prime) = -\frac{[\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]\cdot \vect{n}_{\vect{r}}}{4\pi\sqrt{\det{\bm{\diel}}}
 *     [(\vect{r}-\vect{r}^\prime)^t\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]^{\frac{3}{2}}} 
 *     \end{aligned}\\
 *   &\begin{aligned}
 *     \frac{\partial^2}{\partial{\vect{n}_{\vect{r}}}\partial{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) &=
 *     \frac{\vect{n}_{\vect{r}}\cdot[\bm{\diel}^{-1} \vect{n}_{\vect{r}^\prime}]}{4\pi\sqrt{\det{\bm{\diel}}}
 *     [(\vect{r}-\vect{r}^\prime)^t\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]^{\frac{3}{2}}} \\
 *     &-3\frac{\lbrace[\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]\cdot \vect{n}_{\vect{r}}\rbrace\lbrace[\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]\cdot \vect{n}_{\vect{r}^\prime}\rbrace}{4\pi\sqrt{\det{\bm{\diel}}}
 *     [(\vect{r}-\vect{r}^\prime)^t\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]^{\frac{5}{2}}}
 *    \end{aligned}
 *   \end{align}
 *   \f]
 */
inline Eigen::Array4d analyticAnisotropicLiquid(const Eigen::Vector3d & epsilon, 
    	                    const Eigen::Vector3d & euler,
                                const Eigen::Vector3d & spNormal,
                                const Eigen::Vector3d & sp,
                                const Eigen::Vector3d & ppNormal, const Eigen::Vector3d & pp) {
    Eigen::Array4d result = Eigen::Array4d::Zero();
    Eigen::Matrix3d epsilonInv, R;
    eulerRotation(R, euler);
    Eigen::Vector3d scratch; 
    scratch << (1.0/epsilon(0)), (1.0/epsilon(1)), (1.0/epsilon(2));
    epsilonInv = R * scratch.asDiagonal() * R.transpose();
    double detEps = epsilon(0) * epsilon(1) * epsilon(2);
    Eigen::Vector3d right = epsilonInv * (sp - pp);
    Eigen::Vector3d left  = (sp - pp).transpose();
    double distance   = std::sqrt(left.dot(right));
    double distance_3 = std::pow(distance, 3);
    double distance_5 = std::pow(distance, 5);

    // Value of the function
    result(0) = 1.0 / (std::sqrt(detEps) * distance);
    // Value of the directional derivative wrt probe
    result(1) = right.dot(ppNormal) / (std::sqrt(detEps) * distance_3);
    // Directional derivative wrt source
    result(2) = - right.dot(spNormal) / (std::sqrt(detEps) * distance_3);
    // Value of the Hessian
    Eigen::Vector3d eps_ppNormal = epsilonInv * ppNormal;
    result(3) = spNormal.dot(eps_ppNormal) / (std::sqrt(detEps) * distance_3)
    	  - 3 * (right.dot(spNormal)) * (right.dot(ppNormal)) / (std::sqrt(detEps) * distance_5);

    return result;
}

#endif // ANALYTICEVALUATE_HPP
