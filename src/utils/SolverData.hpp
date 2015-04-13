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

#ifndef SOLVERDATA_HPP
#define SOLVERDATA_HPP

#include "Config.hpp"

class IGreensFunction;

/*! @struct solverData
 *  @brief Contains all data defined from user input in the solver section.
 *  @var solverData::gfInside
 *  The Green's function inside the cavity.
 *  @var solverData::gfOutside
 *  The Green's function outside the cavity.
 *  @var solverData::compressionParameters
 *  A priori and a posteriori compression settings for wavelet solvers.
 *  @var solverData::correction
 *  The correction factor to be use in a CPCM calculation.
 *  @var solverData::integralEquation
 *  The type of integral equation to solve, relevant only for wavelet solvers.
 *  @var solverData::hermitivitize
 *  Triggers hermitivitization of the PCM matrix obtained by collocation.
 */

struct solverData {
    IGreensFunction * gfInside;
    IGreensFunction * gfOutside;
    double correction;
    int integralEquation;
    bool hermitivitize;
    solverData() {}
    solverData(IGreensFunction * _gfInside, IGreensFunction * _gfOutside,
	       double _correction = 0.0,  int _integralEquation = 1, bool _symm = true) :
        gfInside(_gfInside), gfOutside(_gfOutside), correction(_correction),
	integralEquation(_integralEquation), hermitivitize(_symm) {}
};

#endif // SOLVERDATA_HPP
