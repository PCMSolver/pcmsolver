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

#include "GreensFunctionFactory.hpp"

#include <stdexcept>
#include <string>

#include "Config.hpp"

bool GreensFunctionFactory::registerGreensFunction(const std::string & greenID,
        createGreensFunctionCallback createFunction)
{
    return callbacks.insert(CallbackMap::value_type(greenID, createFunction)).second;
}

bool GreensFunctionFactory::unRegisterGreensFunction(const std::string & greenID)
{
    return callbacks.erase(greenID) == 1;
}

IGreensFunction * GreensFunctionFactory::createGreensFunction(
    const std::string & greenID, const greenData & _data)
{
    CallbackMap::const_iterator i = callbacks.find(greenID);
    if (i == callbacks.end()) {
        // The greenID was not found
        throw std::runtime_error("The unknown Green's function ID " + greenID +
                                 " occurred in GreensFunctionFactory.");
    }
    // Invoke the creation function
    return (i->second)(_data);
}
