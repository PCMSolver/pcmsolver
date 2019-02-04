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

#include <cstdio>
#include <ctime>
#include <string>

#include "VersionInfo.hpp"

/*! \file Citation.hpp
 *  \brief Contains the citation text.
 */

inline std::string citation_message() {
  // clang-format off
  const char * fmt =
     "\n-----------------------------------------------------------------------\n"
     "   PCMSolver: An Open Source API for the Polarizable Continuum Model\n"
     "                   PCMSolver %s\n\n"
     "           Git: Branch {%s}, Revision {%s}\n\n"
     " R. Di Remigio, A. H. Steindal, K. Mozgawa, V. Weijo, H. Cao, and\n"
     " L. Frediani, Int. J. Quantum Chem., 2019, 119 (1), e25685.\n\n"
     " Source repository: https://github.com/PCMSolver/pcmsolver\n"
     " Documentation: https://pcmsolver.readthedocs.io/\n"
     " PCMSolver initialized on: %s\n"
     "-----------------------------------------------------------------------\n";
  // clang-format on
  // Get current time
  time_t rawtime;
  struct tm * timeinfo;
  char current_time[80];

  std::time(&rawtime);
  timeinfo = std::localtime(&rawtime);

  std::strftime(
      current_time, sizeof(current_time), "%A, %d %B %Y %I:%M %p", timeinfo);

  char citation[1000];
  std::sprintf(citation,
               fmt,
               PROJECT_VERSION,
               GIT_COMMIT_BRANCH,
               GIT_COMMIT_HASH,
               current_time);

  return std::string(citation);
}
