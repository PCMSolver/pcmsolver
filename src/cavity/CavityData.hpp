/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
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

#include <string>

#include "Config.hpp"

#include "utils/Molecule.hpp"

namespace pcm {
/*! @struct CavityData
 *  @brief Contains all data defined from user input in the cavity section.
 */
struct CavityData {
  /*! The type of cavity */
  std::string cavityType;
  /*! Molecule object with the relevant data for cavity generation */
  Molecule molecule;
  /*! The average tesserae area. Relevant for GePolCavity */
  double area;
  /*! The radius of the spherical probe representing the solvent */
  double probeRadius;
  /*! Triggers the addition of spheres not centered on atoms, relevant for
   * GePolCavity */
  double minimalRadius;
  /*! Name of the .npz file containing the cavity specification for a restart */
  std::string filename;

  CavityData(const std::string & type,
             const Molecule & _molec,
             double _area,
             double _probeRadius,
             double _minRadius,
             const std::string & _fname)
      : cavityType(type),
        molecule(_molec),
        area(_area),
        probeRadius(_probeRadius),
        minimalRadius(_minRadius),
        filename(_fname) {}
};
} // namespace pcm
