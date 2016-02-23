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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#ifndef REGISTERSOLVERTOFACTORY_HPP
#define REGISTERSOLVERTOFACTORY_HPP

#include <string>

#include "SolverData.hpp"
#include "utils/Factory.hpp"
#include "CPCMSolver.hpp"
#include "IEFSolver.hpp"

/*! \file RegisterSolverToFactory.hpp
 *  \brief Register each solver to the factory.
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  This file collects all the creational functions needed for the creation
 *  of the solver objects by means of the factory method.
 *  Originally, each of them was in the same file as the respective class.
 *  This, however, lead to intricate inclusion dependencies.
 */

namespace
{
    PCMSolver * createCPCMSolver(const solverData & data)
    {
        return new CPCMSolver(data.hermitivitize, data.correction);
    }
    const std::string CPCMSOLVER("CPCM");
    const bool registeredCPCMSolver =
        Factory<PCMSolver, solverData>::TheFactory().registerObject(CPCMSOLVER, createCPCMSolver);
}

namespace
{
    PCMSolver * createIEFSolver(const solverData & data)
    {
        return new IEFSolver(data.hermitivitize);
    }
    const std::string IEFSOLVER("IEFPCM");
    const bool registeredIEFSolver =
        Factory<PCMSolver, solverData>::TheFactory().registerObject(IEFSOLVER, createIEFSolver);
}

#endif // REGISTERSOLVERTOFACTORY_HPP
