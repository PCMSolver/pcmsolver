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

/*  To time a code snippet:
 *  \code{.cpp}
 *  TIMER_ON("code-snippet");
 *  // code-snippet
 *  TIMER_OFF("code-snippet");
 *  \endcode
 *  The timings are printed out by a call to the TIMER_DONE function.
 */
#ifdef ENABLE_TIMER

#include "utils/Timer.hpp"

#define TIMER_ON timer::timerON
#define TIMER_OFF timer::timerOFF
#define TIMER_DONE timer::timerDONE

#else /* ENABLE_TIMER */

#define TIMER_ON(...)
#define TIMER_OFF(...)
#define TIMER_DONE(...)

#endif /* ENABLE_TIMER */
