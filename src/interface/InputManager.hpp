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

#ifndef INPUTMANAGER_HPP
#define INPUTMANAGER_HPP

#include <algorithm>
#include <iostream>
#include <string>

#include "Config.hpp"

#include <boost/algorithm/string.hpp>

inline bool invalidChar(char c);

struct cavityInput;
inline void print(const cavityInput &);
inline void init(cavityInput &);
inline void trim(cavityInput &);

struct solverInput;
inline void print(const solverInput &);
inline void init(solverInput &);
inline void trim(solverInput &);

struct greenInput;
inline void print(const greenInput &);
inline void init(greenInput &);
inline void trim(greenInput &);

bool invalidChar(char c)
{
    return !std::isprint( static_cast<unsigned char>( c ) );
}

/*! @struct cavityInput
 *  @brief Data structure for host-API cavity input communication.
 *  @var cavityInput::cavity_type
 *  Type of cavity requested.
 *  @var cavityInput::patch_level
 *  Wavelet cavity mesh patch level.
 *  @var cavityInput::coarsity
 *  Wavelet cavity mesh coarsity.
 *  @var cavityInput::area
 *  Average tesserae area.
 *  @var cavityInput::min_distance
 *  Minimal distance between sampling points.
 *  @var cavityInput::der_order
 *  Derivative order for the switching function.
 *  @var cavityInput::scaling
 *  Whether to scale or not the atomic radii.
 *  @var cavityInput::radii_set
 *  The built-in radii set to be used.
 *  @var cavityInput::restart_name
 *  Name of the .npz file for GePol cavity restart.
 *  @var cavityInput::min_radius
 *  Minimal radius for the added spheres.
 */
struct cavityInput
{
    int patch_level;
    double coarsity;
    double area;
    double min_distance;
    int der_order;
    char cavity_type[8];
    bool scaling;
    char restart_name[20];
    double min_radius;
    char radii_set[8];
};

void print(const cavityInput & inp)
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
}
void init(cavityInput & inp)
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
}
void trim(cavityInput & inp)
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
    */
}

/*! @struct solverInput
 *  @brief Data structure for host-API solver input communication.
 *  @var solverInput::solver_type
 *  Type of solver requested.
 *  @var solverInput::solvent
 *  Name of the solvent.
 *  @var solverInput::equation_type
 *  Type of the integral equation to be used.
 *  @var solverInput::correction
 *  Correction in the CPCM apparent surface charge scaling factor.
 *  @var solverInput::probe_radius
 *  Radius of the spherical probe mimicking the solvent.
 */
struct solverInput
{
    char solver_type[7];
    double correction;
    char solvent[16];
    double probe_radius;
    char equation_type[11];
};

void print(const solverInput & inp)
{
    std::cout << "solver type " << std::string(inp.solver_type) << std::endl;
    std::cout << "solvent " << std::string(inp.solvent) << std::endl;
    std::cout << "equation type " << std::string(inp.equation_type) << std::endl;
    std::cout << "correction " << inp.correction << std::endl;
    std::cout << "probe_radius " << inp.probe_radius << std::endl;
}
void init(solverInput & inp)
{
    std::fill((inp.solver_type), (inp.solver_type) + 7, 0);
    inp.correction = 0.0;
    std::fill((inp.solvent), (inp.solvent) + 16, 0);
    inp.probe_radius = 0.0;
    std::fill((inp.equation_type), (inp.equation_type) + 11, 0);
}
void trim(solverInput & inp)
{
    std::string s1(inp.solver_type);
    s1.erase(std::remove_if(s1.begin(), s1.end(), invalidChar), s1.end());
    boost::algorithm::trim(s1);
    std::fill((inp.solver_type), (inp.solver_type) + 7, 0);
    strncpy(inp.solver_type, s1.c_str(), s1.length());

    std::string s2(inp.solvent);
    s2.erase(std::remove_if(s2.begin(), s2.end(), invalidChar), s2.end());
    boost::algorithm::trim(s2);
    std::fill((inp.solvent), (inp.solvent) + 16, 0);
    strncpy(inp.solvent, s2.c_str(), s2.length());

    std::string s3(inp.equation_type);
    s3.erase(std::remove_if(s3.begin(), s3.end(), invalidChar), s3.end());
    boost::algorithm::trim(s3);
    std::fill((inp.equation_type), (inp.equation_type) + 11, 0);
    strncpy(inp.equation_type, s3.c_str(), s3.length());
}

/*! @struct greenInput
 *  @brief Data structure for host-API Green's function input communication.
 *  @var greenInput::inside_type
 *  Type of Green's function requested inside the cavity.
 *  @var greenInput::outside_epsilon
 *  Value of the static permittivity outside the cavity.
 *  @var greenInput::outside_type
 *  Type of Green's function requested outside the cavity.
 */
struct greenInput
{
    char inside_type[7];
    double outside_epsilon;
    char outside_type[22];
};

void print(const greenInput & inp)
{
    std::cout << "inside type " << std::string(inp.inside_type) << std::endl;
    std::cout << "outside type " << std::string(inp.outside_type) << std::endl;
    std::cout << "epsilon outside " << inp.outside_epsilon << std::endl;
}
void init(greenInput & inp)
{
    std::fill((inp.inside_type), (inp.inside_type) + 7, 0);
    inp.outside_epsilon = 0.0;
    std::fill((inp.outside_type), (inp.outside_type) + 22, 0);
}
void trim(greenInput & inp)
{
    std::string s1(inp.inside_type);
    s1.erase(std::remove_if(s1.begin(), s1.end(), invalidChar), s1.end());
    boost::algorithm::trim(s1);
    std::fill((inp.inside_type), (inp.inside_type) + 7, 0);
    strncpy(inp.inside_type, s1.c_str(), s1.length());

    std::string s2(inp.outside_type);
    s2.erase(std::remove_if(s2.begin(), s2.end(), invalidChar), s2.end());
    boost::algorithm::trim(s2);
    std::fill((inp.outside_type), (inp.outside_type) + 22, 0);
    strncpy(inp.outside_type, s2.c_str(), s2.length());
}

#endif // INPUTMANAGER_HPP
