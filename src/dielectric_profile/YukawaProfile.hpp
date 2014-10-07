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

#ifndef YUKAWAPROFILE_HPP
#define YUKAWAPROFILE_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Dense>

/*! \file YukawaProfile.hpp
 *  \class YukawaProfile
 *  \brief describes a medium with damping, i.e. ionic liquid 
 *  \author Roberto Di Remigio
 *  \date 2014
 */

class YukawaProfile
{
private:
    double epsilon_;
    double kappa_;
public:
    YukawaProfile(double eps, double k) : epsilon_(eps), kappa_(k) {}
    double epsilon() { return epsilon_; }
    double kappa()   { return kappa_; }
    bool isUniform() { return false; }
};

#endif // YUKAWAPROFILE_HPP
