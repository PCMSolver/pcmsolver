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

#include "Element.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <stdexcept>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "Sphere.hpp"

/*! We start from the parametric representation of the sphere:
 * \f[
 * \begin{aligned}
 * (x - x_0) &= R\sin\theta\cos\phi = \sqrt{R^2 - u^2}\cos\phi \\
 * (y - y_0) &= R\sin\theta\sin\phi = \sqrt{R^2 - u^2}\sin\phi \\
 * (z - z_0) &= R\cos\theta         = u
 * \end{aligned}
 * \f]
 * \f$\phi\f$ is called azimuth or azimuthal angle, while \f$\theta\f$ is the polar angle
 * (physicists' convention)
 *
 * Their inverse:
 * \f[
 * \begin{aligned}
 * R      &= \sqrt{(x-x_0)^2 + (y-y_0)^2 + (z-z_0)^2} \\
 * \phi   &= \arctan{\frac{(y-y_0)}{(x-x_0)}}         \\
 * \theta &= \arccos{\frac{(z-z_0)}{R}}
 * \end{aligned}
 * \f]
 * Given the parametric description above, the following vectors
 * span the tangent space to the sphere (might be that they are undefined
 * in some points...)
 * \f[
 * v_\phi = \begin{pmatrix}
 *        \partial_\theta x \\
 *        \partial_\theta y \\
 *        \partial_\theta z
 *       \end{pmatrix}
 *       = \begin{pmatrix}
 *         -\sqrt{R^2-u^2}\sin\phi \\
 *          \sqrt{R^2-u^2}\cos\phi \\
 *          0
 *       \end{pmatrix}
 * \quad\quad
 * v_u = \begin{pmatrix}
 *        \partial_u x \\
 *        \partial_u y \\
 *        \partial_u z 
 *       \end{pmatrix}
 *     = \begin{pmatrix}
 *         -\frac{u}{\sqrt{R^2-u^2}}\cos\phi \\
 *         -\frac{u}{\sqrt{R^2-u^2}}\sin\phi \\
 *          1
 *       \end{pmatrix}
 *  \f]
 *
 * The tangent and bitangent are then:
 * \f[
 * \begin{alignat}{2}
 * t = -\frac{v_\phi}{|v_\phi|}, \quad&\quad
 * b = -\frac{v_u}{|v_u|}
 * \end{alignat}
 * \f]
 * the minus sign is required to preserve consistency with GAMESS subroutines
 */
void Element::tangent_and_bitangent(Eigen::Vector3d & t_, Eigen::Vector3d & b_) const
{
	double polar_angle = std::acos((center_(2) - sphere_.center(2)) / sphere_.radius()); // \theta
	double azimuth = std::atan((center_(1) - sphere_.center(1)) / (center_(0) - sphere_.center(0))); // \phi
	double Rcos = sphere_.radius() * std::cos(polar_angle); // R\cos\theta
	double Rsin = std::sqrt(std::pow(sphere_.radius(), 2) - std::pow(Rcos, 2)); // R\sin\theta
	t_ << Rsin*std::sin(azimuth), -Rsin*std::cos(azimuth), 0.0;
	t_.normalize();
	b_ << (Rcos/Rsin)*std::cos(azimuth), (Rcos/Rsin)*std::sin(azimuth), -1.0;
	b_.normalize();
}
