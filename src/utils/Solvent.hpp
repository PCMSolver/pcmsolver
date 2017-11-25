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

#pragma once

#include <iosfwd>
#include <map>
#include <string>

#include "Config.hpp"

/*! \file Solvent.hpp */

namespace pcm {
namespace utils {

/*! \struct Solvent
 *  \brief POD describing a solvent.
 *  \author Roberto Di Remigio
 *  \date 2011, 2016
 *
 * A Solvent object contains all the solvent-related experimental data
 * needed to set up the Green's functions and the non-electrostatic
 * terms calculations.
 */

struct Solvent {
  Solvent() {}
  Solvent(const std::string & n, double es, double ed, double r)
      : name(n), epsStatic(es), epsDynamic(ed), probeRadius(r) {}
  /*! Solvent name */
  std::string name;
  /*! Static permittivity, in AU*/
  double epsStatic;
  /*! Optical permittivity, in AU */
  double epsDynamic;
  /*! Radius of the spherical probe mimicking the solvent, in Angstrom */
  double probeRadius;
};
} // namespace utils

std::ostream & operator<<(std::ostream & os, utils::Solvent & solvent);

/*! \brief typedef for the map between solvent name and Solvent object. */
typedef std::map<std::string, utils::Solvent> SolventMap;

/*! \brief Returns the map between solvent names and Solvent objects.
 *
 *  This map contains solvent data taken from the DALTON2011 internal
 *  implementation of the Polarizable Continuum Model.
 */
SolventMap & solvents();
} // namespace pcm
