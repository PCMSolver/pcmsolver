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

#include "DiagonalIntegratorFactory.hpp"

#include <stdexcept>
#include <string>

#include "Config.hpp"

bool DiagonalIntegratorFactory::registerDiagonalIntegrator(const std::string & integratorID,
        createDiagonalIntegratorCallback createFunction)
{
    return callbacks.insert(CallbackMap::value_type(integratorID, createFunction)).second;
}

bool DiagonalIntegratorFactory::unRegisterDiagonalIntegrator(const std::string & integratorID)
{
    return callbacks.erase(integratorID) == 1;
}

DiagonalIntegrator * DiagonalIntegratorFactory::createDiagonalIntegrator(
    const std::string & integratorID)
{
    CallbackMap::const_iterator i = callbacks.find(integratorID);
    if (i == callbacks.end()) {
        // The integratorID was not found
        throw std::runtime_error("The unknown Green's function ID " + integratorID +
                                 " occurred in DiagonalIntegratorFactory.");
    }
    // Invoke the creation function
    return (i->second)();
}
