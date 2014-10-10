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

#ifndef CITATION_HPP
#define CITATION_HPP

#include <sstream>
#include <string>

#include "Config.hpp"

// This is to stringify the PROJECT_VERSION preprocessor constant
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)

/*! \file Citation.hpp
 *  \brief Contains the citation text.
 *
 */

inline std::string citation_message()
{
    std::ostringstream rest;
    std::string version(TOSTRING(PROJECT_VERSION));
    rest << "\n" << std::endl;
    rest << " * PCMSolver, an API for the Polarizable Continuum Model electrostatic problem. Version "
         << version << std::endl;
    rest << "   Main authors: R. Di Remigio, L. Frediani, K. Mozgawa" << std::endl;
    rest << "    With contributions from:" << std::endl;
    rest << "     R. Bast            (CMake framework)" << std::endl;
#if defined (DEVELOPMENT_CODE)
    rest << "     S. Corni           (Metal sphere Green's function)" << std::endl;
#endif
    rest << "     U. Ekstroem        (automatic differentiation library)" << std::endl;
#if defined (WAVELET_DEVELOPMENT)
    rest << "     H. Harbrecht       (wavelet cavity and solvers libraries)" <<
         std::endl;
#endif
    rest << "     J. Juselius        (input parsing library and CMake framework)" <<
         std::endl;
#if defined (TSLESS_DEVELOPMENT)
    rest << "     C. S. Pomelli      (TsLess cavity library)" << std::endl;
#endif
#if defined (WAVELET_DEVELOPMENT)
    rest << "     M. Randrianarivony (wavelet cavity library)" << std::endl;
    rest << "     V. Weijo           (wavelet libraries and cavity visualization scripts)"
         << std::endl;
#endif
    rest << "   Theory: - J. Tomasi, B. Mennucci and R. Cammi:" << std::endl;
    rest << "            \"Quantum Mechanical Continuum Solvation Models\", Chem. Rev., 105 (2005) 2999"
         << std::endl;
    rest << "   PCMSolver is distributed under the terms of the GNU Lesser General Public License."
         << std::endl;
    return rest.str();
}
#endif // CITATION_HPP
