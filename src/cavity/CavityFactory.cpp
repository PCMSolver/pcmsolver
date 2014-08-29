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
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#include "CavityFactory.hpp"

#include <string>
#include <stdexcept>
#include <vector>

#include "Config.hpp"


#include "Sphere.hpp"

bool CavityFactory::registerCavity(std::string cavityID,
                                   createCavityCallback createFunction)
{
    return callbacks.insert(CallbackMap::value_type(cavityID, createFunction)).second;
}

bool CavityFactory::unRegisterCavity(std::string cavityID)
{
    return callbacks.erase(cavityID) == 1;
}

Cavity * CavityFactory::createCavity(std::string cavityID, const cavityData & _data)
{
    CallbackMap::const_iterator i = callbacks.find(cavityID);
    if (i == callbacks.end()) {
        // The cavityID was not found
        throw std::runtime_error("The unknown cavity ID " + cavityID +
                                 " occurred in CavityFactory.");
    }
    // Invoke the creation function
    return (i->second)(_data);
}
