/**
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

#pragma once

#include <cmath>

/*! \file PhysicalConstants.hpp
 *  \brief This header file contains physical constants to be used throughout the
 * module.
 *  \author Roberto Di Remigio
 *  \date 2013
 */

struct LengthConversion {
  double operator()() { return BOHR_TO_ANGSTROM; }
  double BOHR_TO_ANGSTROM;
};

extern LengthConversion bohrToAngstrom;

void initBohrToAngstrom(LengthConversion & conversion, int year = 2010);

double angstromToBohr();

double bohr2ToAngstrom2();

double angstrom2ToBohr2();

double bohr3ToAngstrom3();

double angstrom3ToBohr3();
