/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#include "PhysicalConstants.hpp"

LengthConversion bohrToAngstrom;

void initBohrToAngstrom(LengthConversion & conversion, int year) {
  switch (year) {
    case 2010:
      /*
       * P. J. Mohr, B. N. Taylor, and D. B. Newell, Rev. Mod. Phys. 84(4), 1527-1605
       * (2012)
       */
      conversion.BOHR_TO_ANGSTROM = 0.52917721092;
      break;
    case 2006:
      /*
       * P. J. Mohr, B. N. Taylor, and D. B. Newell, Rev. Mod. Phys. 80(2), 633-730
       * (2008)
       */
      conversion.BOHR_TO_ANGSTROM = 0.52917720859;
      break;
    case 2002:
      /* Choose this to match DALTON values for physical constants
       * P. J. Mohr and B. N. Taylor, Rev. Mod. Phys. 77(1), 1-107 (2005)
       */
      conversion.BOHR_TO_ANGSTROM = 0.5291772108;
      break;
    case 1998:
      /*
       * P. J. Mohr and B. N. Taylor, Rev. Mod. Phys. 72(2), 351-495 (2000)
       */
      conversion.BOHR_TO_ANGSTROM = 0.5291772083;
      break;
    default:
      /*
       * Use latest values by default */
      conversion.BOHR_TO_ANGSTROM = 0.52917721092;
      break;
  }
}

double angstromToBohr() { return (1.0 / bohrToAngstrom()); }

double bohr2ToAngstrom2() { return std::pow(bohrToAngstrom(), 2); }

double angstrom2ToBohr2() { return (1.0 / bohr2ToAngstrom2()); }

double bohr3ToAngstrom3() { return std::pow(bohrToAngstrom(), 3); }

double angstrom3ToBohr3() { return (1.0 / bohr3ToAngstrom3()); }
