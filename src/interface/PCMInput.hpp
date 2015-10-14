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

#include <algorithm>
#include <iostream>
#include <string>

#include "Config.hpp"

#include <boost/algorithm/string.hpp>

inline bool invalidChar(char c);

struct PCMInput;
inline void print(const PCMInput &);
inline void init(PCMInput &);
inline void trim(PCMInput &);

bool invalidChar(char c)
{
    return !std::isprint( static_cast<unsigned char>( c ) );
}

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
	/// Minimal distance between sampling points.
	double min_distance;
	/// Derivative order for the switching function.
	int der_order;
	/// Whether to scale or not the atomic radii.
	bool scaling;
	/// The built-in radii set to be used.
	char radii_set[8];
	/// Name of the .npz file for GePol cavity restart.
	char restart_name[20];
	/// Minimal radius for the added spheres.
	double min_radius;
	/// Type of solver requested.
	char solver_type[7];
	/// Name of the solvent.
	char solvent[16];
	/// Type of the integral equation to be used.
	char equation_type[11];
	/// Correction in the CPCM apparent surface charge scaling factor.
	double correction;
	/// Radius of the spherical probe mimicking the solvent.
	double probe_radius;
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
void init(PCMInput & inp)
{
    inp.patch_level = 0;
    std::fill((inp.cavity_type), (inp.cavity_type) + 8, 0);
    inp.coarsity = 0.0;
    inp.area = 0.0;
    inp.min_distance = 0.0;
    inp.der_order = 0;
    inp.scaling = false;
    std::fill((inp.radii_set), (inp.radii_set) + 8, 0);
    std::fill((inp.restart_name), (inp.restart_name) + 20, 0);
    inp.min_radius = 0.0;
    std::fill((inp.solver_type), (inp.solver_type) + 7, 0);
    inp.correction = 0.0;
    std::fill((inp.solvent), (inp.solvent) + 16, 0);
    inp.probe_radius = 0.0;
    std::fill((inp.equation_type), (inp.equation_type) + 11, 0);
    std::fill((inp.inside_type), (inp.inside_type) + 7, 0);
    inp.outside_epsilon = 0.0;
    std::fill((inp.outside_type), (inp.outside_type) + 22, 0);
}
void trim(PCMInput & inp)
{
    std::string s1(inp.cavity_type);
    s1.erase(std::remove_if(s1.begin(), s1.end(), invalidChar), s1.end());
    boost::algorithm::trim(s1);
    std::fill((inp.cavity_type), (inp.cavity_type) + 8, 0);
    strncpy(inp.cavity_type, s1.c_str(), s1.length());

    /*
    std::string s2(inp.restart_name);
    s2.erase(std::remove_if(s2.begin(), s2.end(), invalidChar), s2.end());
    boost::algorithm::trim(s2);
    std::fill((inp.restart_name), (inp.restart_name) + 20, 0);
    strncpy(inp.cavity_type, s2.c_str(), s2.length());

    std::string s3(inp.radii_set);
    s3.erase(std::remove_if(s3.begin(), s3.end(), invalidChar), s3.end());
    boost::algorithm::trim(s3);
    std::fill((inp.radii_set), (inp.radii_set) + 8, 0);
    strncpy(inp.cavity_type, s3.c_str(), s3.length());

    std::string s4(inp.solver_type);
    s4.erase(std::remove_if(s4.begin(), s4.end(), invalidChar), s4.end());
    boost::algorithm::trim(s4);
    std::fill((inp.solver_type), (inp.solver_type) + 7, 0);
    strncpy(inp.solver_type, s4.c_str(), s4.length());

    std::string s5(inp.solvent);
    s5.erase(std::remove_if(s5.begin(), s5.end(), invalidChar), s5.end());
    boost::algorithm::trim(s5);
    std::fill((inp.solvent), (inp.solvent) + 16, 0);
    strncpy(inp.solvent, s5.c_str(), s5.length());

    std::string s6(inp.equation_type);
    s6.erase(std::remove_if(s6.begin(), s6.end(), invalidChar), s6.end());
    boost::algorithm::trim(s6);
    std::fill((inp.equation_type), (inp.equation_type) + 11, 0);
    strncpy(inp.equation_type, s6.c_str(), s6.length());

    std::string s7(inp.inside_type);
    s7.erase(std::remove_if(s7.begin(), s7.end(), invalidChar), s7.end());
    boost::algorithm::trim(s7);
    std::fill((inp.inside_type), (inp.inside_type) + 7, 0);
    strncpy(inp.inside_type, s7.c_str(), s7.length());

    std::string s8(inp.outside_type);
    s8.erase(std::remove_if(s8.begin(), s8.end(), invalidChar), s8.end());
    boost::algorithm::trim(s8);
    std::fill((inp.outside_type), (inp.outside_type) + 22, 0);
    strncpy(inp.outside_type, s8.c_str(), s8.length());
    */
}

#endif // PCMINPUT_HPP
