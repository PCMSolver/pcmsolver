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
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#ifndef UNIFORMPROFILE_HPP
#define UNIFORMPROFILE_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Dense>

/*! \file UniformProfile.hpp
 *  \class UniformProfile
 *  \brief describes a uniform dielectric profile
 *  \author Roberto Di Remigio
 *  \date 2014
 */

class UniformProfile
{
private:
    double epsilon_;
public:
    UniformProfile(double eps) : epsilon_(eps) {}
    double epsilon() { return epsilon_; }
    bool isUniform() { return true; }
};

#endif // UNIFORMPROFILE_HPP
