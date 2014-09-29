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

#ifndef INPUTMANAGER_HPP
#define INPUTMANAGER_HPP

#include <iostream>
#include <string>

#include "Config.hpp"

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
 *  @var min_distance
 *  Minimal distance between sampling points.
 *  @var der_order
 *  Derivative order for the switching function.
 *  @var scaling
 *  Whether to scale or not the atomic radii.
 *  @var radii_set
 *  The built-in radii set to be used.
 *  @var restart_name
 *  Name of the .npz file for GePol cavity restart.
 *  @var min_radius
 *  Minimal radius for the added spheres.
 */
struct cavityInput {
    char cavity_type[8];
    int patch_level;
    double coarsity;
    double area;
    double min_distance;
    int der_order;
    bool scaling;
    char radii_set[8];
    int pad; // Without padding the two strings are lumped together...
    char restart_name[20];
    double min_radius;
    friend std::ostream & operator<<(std::ostream & os, cavityInput & o) {
        os << "cavity type " << std::string(o.cavity_type) << std::endl;
        os << "patch level " << o.patch_level << std::endl;
        os << "coarsity " << o.coarsity << std::endl;
        os << "area " << o.area << std::endl;
        os << "min distance " << o.min_distance << std::endl;
        os << "der order " << o.der_order << std::endl;
        os << "scaling " << o.scaling << std::endl;
        os << "radii set " << std::string(o.radii_set) << std::endl;
        os << "restart name " << std::string(o.restart_name) << std::endl;
	os << "min radius " << o.min_radius;
        return os;
    }
};

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
struct solverInput {
    char solver_type[7];
    int pad1; // Without padding the two strings are lumped together...
    char solvent[16];
    char equation_type[11];
    double correction;
    double probe_radius;
    friend std::ostream & operator<<(std::ostream & os, solverInput & o) {
        os << "solver type " << std::string(o.solver_type) << std::endl;
        os << "solvent " << std::string(o.solvent) << std::endl;
	os << "equation type " << std::string(o.equation_type) << std::endl;
	os << "correction " << o.correction << std::endl;
	os << "probe_radius " << o.probe_radius;
        return os;
    }
};

/*! @struct greenInput
 *  @brief Data structure for host-API Green's function input communication.
 *  @var greenInput::inside_type
 *  Type of Green's function requested inside the cavity.
 *  @var greenInput::outside_epsilon
 *  Value of the static permittivity outside the cavity.
 *  @var greenInput::outside_type
 *  Type of Green's function requested outside the cavity.
 */
struct greenInput {
    char inside_type[7];
    double outside_epsilon;
    char outside_type[22];
    friend std::ostream & operator<<(std::ostream & os, greenInput & o) {
        os << "inside type " << std::string(o.inside_type) << std::endl;
        os << "outside type " << std::string(o.outside_type) << std::endl;
	os << "epsilon outside " << o.outside_epsilon;
        return os;
    }
};

#endif // INPUTMANAGER_HPP
