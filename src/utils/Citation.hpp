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
    rest << " * PCMSolver, an API for the Polarizable Continuum Model electrostatic problem. Version " << version << "\n";
    rest << "   Main authors: R. Di Remigio, L. Frediani, K. Mozgawa\n";
    rest << "    With contributions from:\n";
    rest << "     R. Bast            (CMake framework)\n";                           
    rest << "     U. Ekstroem        (automatic differentiation library)\n";
    rest << "     H. Harbrecht       (wavelet cavity and solvers libraries)\n";
    rest << "     J. Juselius        (input parsing library and CMake framework)\n";
    rest << "     M. Randrianarivony (wavelet cavity library)\n";
    rest << "     V. Weijo           (wavelet libraries and cavity visualization scripts)\n";
    rest << "   Theory: - J. Tomasi, B. Mennucci and R. Cammi:\n";
    rest << "            \"Quantum Mechanical Continuum Solvation Models\", Chem. Rev., 105 (2005) 2999\n";
    rest << "   PCMSolver is distributed under the terms of the GNU Lesser General Public License.\n";
    return rest.str();
}
#endif // CITATION_HPP
