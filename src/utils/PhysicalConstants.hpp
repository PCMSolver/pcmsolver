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

#ifndef PHYSICALCONSTANTS_HPP
#define PHYSICALCONSTANTS_HPP

#include <cmath>

#include "Config.hpp"

/*! \file PhysicalConstants.hpp
 *  \brief This header file contains physical constants to be used throughout the module.
 *  \author Roberto Di Remigio
 *  \date 2013
 */

inline double bohrToAngstrom(int year)
{
    double bohrToAngstrom_;
    switch(year) {
    case 2010:
        /*
         * P. J. Mohr, B. N. Taylor, and D. B. Newell, Rev. Mod. Phys. 84(4), 1527-1605 (2012)
        */
        bohrToAngstrom_ = 0.52917721092;
        break;
    case 2006:
        /*
         * P. J. Mohr, B. N. Taylor, and D. B. Newell, Rev. Mod. Phys. 80(2), 633-730 (2008)
         */
        bohrToAngstrom_ = 0.52917720859;
        break;
    case 2002:
        /* Choose this to match DALTON values for physical constants
         * P. J. Mohr and B. N. Taylor, Rev. Mod. Phys. 77(1), 1-107 (2005)
         */
        bohrToAngstrom_ = 0.5291772108;
        break;
    case 1998:
        /*
         * P. J. Mohr and B. N. Taylor, Rev. Mod. Phys. 72(2), 351-495 (2000)
         */
        bohrToAngstrom_ = 0.5291772083;
        break;
    default:
        /*
         * Use latest values by default
         */
        bohrToAngstrom_ = 0.52917721092;
        break;
    }

    return bohrToAngstrom_;
}

inline double angstromToBohr(int year)
{
    return (1.0 / bohrToAngstrom(year));
}

inline double bohr2ToAngstrom2(int year)
{
    return std::pow(bohrToAngstrom(year), 2);
}

inline double angstrom2ToBohr2(int year)
{
    return (1.0 / bohr2ToAngstrom2(year));
}

inline double bohr3ToAngstrom3(int year)
{
    return std::pow(bohrToAngstrom(year), 3);
}

inline double angstrom3ToBohr3(int year)
{
    return (1.0 / bohr3ToAngstrom3(year));
}

/*
 *  Taken from the CODATA 2010 recommended values set:
 *  http://physics.nist.gov/cuu/Constants/index.html
 *  Left for backward-compatibility in the unit tests!
 */
/*! \brief Bohr to Angstrom conversion factor
 */
static const double convertBohrToAngstrom = 0.52917721092;
/*! \brief Bohr^2 to Angstrom^2 conversion factor
 */
static const double convertBohr2ToAngstrom2 = 0.2800285205570702;
/*! \brief Bohr^3 to Angstrom^3 conversion factor
 */
static const double convertBohr3ToAngstrom3 = 0.14818471148644433;

#endif // PHYSICALCONSTANTS_HPP
