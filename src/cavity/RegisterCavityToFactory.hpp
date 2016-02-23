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

#ifndef REGISTERCAVITYTOFACTORY_HPP
#define REGISTERCAVITYTOFACTORY_HPP

#include <string>

#include "CavityData.hpp"
#include "utils/Factory.hpp"
#include "GePolCavity.hpp"
#include "RestartCavity.hpp"

/*! \file RegisterCavityToFactory.hpp
 *  \brief Register each cavity to the factory.
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  This file collects all the creational functions needed for the creation
 *  of the cavity objects by means of the factory method.
 *  Originally, each of them was in the same file as the respective class.
 *  This, however, lead to intricate inclusion dependencies.
 */

namespace
{
    Cavity * createGePolCavity(const cavityData & data)
    {
        return new GePolCavity(data.molecule, data.area, data.probeRadius,
                               data.minimalRadius);
    }
    const std::string GEPOL("GEPOL");
    const bool registeredGePol = Factory<Cavity, cavityData>::TheFactory().registerObject(GEPOL,
                                 createGePolCavity);
}

namespace
{
    Cavity * createRestartCavity(const cavityData & data)
    {
        return new RestartCavity(data.filename);
    }
    const std::string RESTART("RESTART");
    const bool registeredRestart = Factory<Cavity, cavityData>::TheFactory().registerObject(
                                       RESTART, createRestartCavity);
}

#endif // REGISTERCAVITYTOFACTORY_HPP
