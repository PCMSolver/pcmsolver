/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#ifndef UNIFORM_HPP
#define UNIFORM_HPP

#include <iostream>

#include "Config.hpp"


/*! \file Uniform.hpp
 *  \struct Uniform
 *  \brief a uniform dielectric profile
 *  \author Roberto Di Remigio
 *  \date 2015
 */

struct Uniform __final
{
    double epsilon;
    Uniform() : epsilon(1.0) {}
    Uniform(double e) : epsilon(e) {}
    friend std::ostream & operator<<(std::ostream & os, Uniform & arg) {
        os << "Permittivity = " << arg.epsilon;
        return os;
    }
};

#endif // UNIFORM_HPP
