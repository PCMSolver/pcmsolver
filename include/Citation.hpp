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

#include <algorithm>
#include <sstream>
#include <string>

#include "Config.hpp"
#include "GitInfo.hpp"

// This is to stringify the PROJECT_VERSION preprocessor constant
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)

/*! \file Citation.hpp
 *  \brief Contains the citation text.
 */

inline std::string citation_message() {
  std::ostringstream rest;
  std::string version(TOSTRING(PROJECT_VERSION));
  rest << "\n" << std::endl;
  rest << " * PCMSolver, an API for the Polarizable Continuum Model electrostatic "
          "problem. Version "
       << version << std::endl;
  rest << "   Main authors: R. Di Remigio, L. Frediani, K. Mozgawa" << std::endl;
  rest << "    With contributions from:" << std::endl;
  rest << "     R. Bast            (CMake framework)" << std::endl;
  rest << "     U. Ekstroem        (automatic differentiation library)" << std::endl;
  rest << "     J. Juselius        (input parsing library and CMake framework)"
       << std::endl;
  rest << "   Theory: - J. Tomasi, B. Mennucci and R. Cammi:" << std::endl;
  rest << "            \"Quantum Mechanical Continuum Solvation Models\", Chem. "
          "Rev., 105 (2005) 2999"
       << std::endl;
  rest << "   PCMSolver is distributed under the terms of the GNU Lesser General "
          "Public License."
       << std::endl;
  return rest.str();
}

inline std::string version_info() {
  std::ostringstream retval;
  retval << " * Git last commit hash   : " << GIT_COMMIT_HASH << std::endl;
  retval << " * Git last commit date   : " << GIT_COMMIT_DATE << std::endl;
  retval << " * Git last commit author : " << GIT_COMMIT_AUTHOR << std::endl;
  return retval.str();
}
