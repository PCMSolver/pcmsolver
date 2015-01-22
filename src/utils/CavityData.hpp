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

#ifndef CAVITYDATA_HPP
#define CAVITYDATA_HPP

#include <string>

#include "Config.hpp"

#include "Molecule.hpp"

/*! @struct cavityData
 *  @brief Contains all data defined from user input in the cavity section.
 *  @var cavityData::molecule
 *  Molecule object with the relevant data for cavity generation.
 *  @var cavityData::area
 *  The average tesserae area. Relevant for GePolCavity.
 *  @var cavityData::probeRadius
 *  The radius of the spherical probe representing the solvent.
 *  @var cavityData::minDistance
 *  The minimal distance between two sampling
 *  points on different spheres. Relevant for TsLessCavity.
 *  @var cavityData::derOrder
 *  The maximum derivative order to be used in the definition
 *  of the smoothing function. Relevant for TsLessCavity.
 *  @var cavityData::minimalRadius
 *  Triggers the addition of spheres not centered on atoms.
 *  Relevant for GePolCavity.
 *  @var cavityData::patchLevel
 *  Relevant for WaveletCavity.
 *  @var cavityData::coarsity
 *  Relevant for WaveletCavity.
 *  @var cavityData::filename
 *  Name of the file containing the cavity
 *  specification for a restart.
 *  @var cavityData::dyadicFile
 */

struct cavityData {
    Molecule molecule;
    double area;
    double probeRadius;
    double minDistance;
    int derOrder;
    double minimalRadius;
    int patchLevel;
    double coarsity;
    std::string filename;
    std::string dyadicFile;
    bool empty;
    cavityData() { empty = true; }
    cavityData(const Molecule & _molec, double _area, double _probeRadius,
               double _minDistance, int _derOrder, double _minRadius,
               int _patchLevel, double _coarsity, const std::string & _fname,
	       const std::string & _dyad) :
        molecule(_molec), area(_area), probeRadius(_probeRadius),
        minDistance(_minDistance), derOrder(_derOrder), minimalRadius(_minRadius),
        patchLevel(_patchLevel), coarsity(_coarsity), filename(_fname),
	dyadicFile(_dyad) { empty = false; }
};

#endif // CAVITYDATA_HPP
