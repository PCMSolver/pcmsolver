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

#include <boost/algorithm/string.hpp>

/*! @struct cavityInput
 *  @brief Data structure for host-API cavity input communication. 
 *  @var cavityInput::type
 *  Type of cavity requested. 
 *  @var cavityInput::area
 *  Average tesserae area. 
 */
struct cavityInput 
{
	char cavity_type[8];
	int patch_level;
	double coarsity;
	double area;
	double min_distance;
	int der_order;
	bool scaling;
	char radii_set[8];
	char restart_name[20];
    friend std::ostream & operator<<(std::ostream & os, cavityInput & o) {
        os << "cavity type " << std::string(o.cavity_type) << std::endl;
	os << "patch level " << o.patch_level << std::endl;
	os << "coarsity " << o.coarsity << std::endl;
	os << "area " << o.area << std::endl;
	os << "min distance " << o.min_distance << std::endl;
	os << "der order " << o.der_order << std::endl;
	os << "scaling " << o.scaling << std::endl;
	os << "radii set " << std::string(o.radii_set) << std::endl;
	os << "restart name " << std::string(o.restart_name);
        return os;
    }
};

/*! @struct solverInput
 *  @brief Data structure for host-API solver input communication. 
 *  @var solverInput::type
 *  Type of solver requested. 
 */
struct solverInput 
{
	char solver_type[7];
	char solvent[16];
    friend std::ostream & operator<<(std::ostream & os, solverInput & o) {
        os << "solver type " << std::string(o.solver_type) << std::endl;
	os << "solvent " << std::string(o.solvent);
        return os;
    }
};

/*! @struct greenInput
 *  @brief Data structure for host-API Green's function input communication. 
 *  @var greenInput::inside_type
 *  Type of Green's function requested inside the cavity.
 *  @var greenInput::outside_type
 *  Type of Green's function requested outside the cavity.
 */
struct greenInput 
{
	int inside_type;
	int outside_type;
    friend std::ostream & operator<<(std::ostream & os, greenInput & o) {
        os << "inside type " << o.inside_type << std::endl;
        os << "outside type " << o.outside_type;
        return os;
    }
};

#endif // INPUTMANAGER_HPP
