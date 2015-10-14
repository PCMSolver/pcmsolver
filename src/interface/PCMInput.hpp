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
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef PCMINPUT_HPP
#define PCMINPUT_HPP

#include <iostream>
#include <string>

#include "Config.hpp"

struct PCMInput;
inline void print(const PCMInput &);

/*! @struct PCMInput
 *  @brief Data structure for host-API input communication.
 */
struct PCMInput
{
	/// Type of cavity requested.
	char cavity_type[8];
	/// Wavelet cavity mesh patch level.
	int patch_level;
	/// Wavelet cavity mesh coarsity.
	double coarsity;
	/// Average tesserae area.
	double area;
	/// The built-in radii set to be used.
	char radii_set[8];
	/// Minimal distance between sampling points.
	double min_distance;
	/// Derivative order for the switching function.
	int der_order;
	/// Whether to scale or not the atomic radii.
	bool scaling;
	/// Name of the .npz file for GePol cavity restart.
	char restart_name[20];
	/// Minimal radius for the added spheres.
	double min_radius;
	/// Type of solver requested.
	char solver_type[7];
	/// Correction in the CPCM apparent surface charge scaling factor.
	double correction;
	/// Name of the solvent.
	char solvent[16];
	/// Radius of the spherical probe mimicking the solvent.
	double probe_radius;
	/// Type of the integral equation to be used.
	char equation_type[11];
	/// Type of Green's function requested inside the cavity.
	char inside_type[7];
	/// Value of the static permittivity outside the cavity.
	double outside_epsilon;
	/// Type of Green's function requested outside the cavity.
	char outside_type[22];
};

void print(const PCMInput & inp)
{
    std::cout << "cavity type " << std::string(inp.cavity_type) << std::endl;
    std::cout << "patch level " << inp.patch_level << std::endl;
    std::cout << "coarsity " << inp.coarsity << std::endl;
    std::cout << "area " << inp.area << std::endl;
    std::cout << "min distance " << inp.min_distance << std::endl;
    std::cout << "der order " << inp.der_order << std::endl;
    std::cout << "scaling " << inp.scaling << std::endl;
    std::cout << "radii set " << std::string(inp.radii_set) << std::endl;
    std::cout << "restart name " << std::string(inp.restart_name) << std::endl;
    std::cout << "min radius " << inp.min_radius << std::endl;
    std::cout << "solver type " << std::string(inp.solver_type) << std::endl;
    std::cout << "solvent " << std::string(inp.solvent) << std::endl;
    std::cout << "equation type " << std::string(inp.equation_type) << std::endl;
    std::cout << "correction " << inp.correction << std::endl;
    std::cout << "probe_radius " << inp.probe_radius << std::endl;
    std::cout << "inside type " << std::string(inp.inside_type) << std::endl;
    std::cout << "outside type " << std::string(inp.outside_type) << std::endl;
    std::cout << "epsilon outside " << inp.outside_epsilon << std::endl;
}

#endif // PCMINPUT_HPP
