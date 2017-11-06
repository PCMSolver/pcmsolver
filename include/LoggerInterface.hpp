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

#ifdef ENABLE_LOGGER

#include "utils/Logger.hpp"
#include "utils/Timer.hpp"

static logging::logger<logging::FileLogPolicy> loggerInstance(
    "pcmsolver.execution.log");

#define LOG loggerInstance.print<logging::printLevel::coarse>
#define LOG_FINE loggerInstance.print<logging::printLevel::fine>
#define LOG_ALL loggerInstance.print<logging::printLevel::everything>
#define LOG_TIME                                                                    \
  loggerInstance.print<logging::printLevel::timings>(timer::Timer::TheTimer());

#else /* ENABLE_LOGGER */

#define LOG(...)
#define LOG_FINE(...)
#define LOG_ALL(...)
#define LOG_TIME

#endif /* ENABLE_LOGGER */
