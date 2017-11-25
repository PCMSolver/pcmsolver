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

#include <string>

#include "Config.hpp"

#include "utils/Molecule.hpp"

/*! @struct CavityData
 *  @brief Contains all data defined from user input in the cavity section.
 */

namespace pcm {
struct CavityData {
  /*! Molecule object with the relevant data for cavity generation */
  Molecule molecule;
  /*! The average tesserae area. Relevant for GePolCavity */
  double area;
  /*! The radius of the spherical probe representing the solvent */
  double probeRadius;
  /*! The minimal distance between two sampling points on different spheres. Relevant
   * for TsLessCavity */
  double minDistance;
  /*! The maximum derivative order to be used in the definition of the smoothing
   * function. Relevant for TsLessCavity */
  int derOrder;
  /*! Triggers the addition of spheres not centered on atoms, relevant for
   * GePolCavity */
  double minimalRadius;
  /*! Patch level, relevant for WaveletCavity */
  int patchLevel;
  /*! Relevant for WaveletCavity */
  double coarsity;
  /*! Name of the .npz file containing the cavity specification for a restart */
  std::string filename;
  /*! Name of the text file with list of points for a wavelet cavityr restart */
  std::string dyadicFile;
  /*! Whether the structure was initialized with user input or not */
  bool empty;

  CavityData() { empty = true; }
  CavityData(const Molecule & _molec,
             double _area,
             double _probeRadius,
             double _minDistance,
             int _derOrder,
             double _minRadius,
             int _patchLevel,
             double _coarsity,
             const std::string & _fname,
             const std::string & _dyad)
      : molecule(_molec),
        area(_area),
        probeRadius(_probeRadius),
        minDistance(_minDistance),
        derOrder(_derOrder),
        minimalRadius(_minRadius),
        patchLevel(_patchLevel),
        coarsity(_coarsity),
        filename(_fname),
        dyadicFile(_dyad) {
    empty = false;
  }
};
} // namespace pcm
