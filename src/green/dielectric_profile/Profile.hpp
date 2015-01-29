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

#ifndef PROFILE_HPP
#define PROFILE_HPP

#include <iosfwd>

#include "Config.hpp"

/*! \file Profile.hpp
 *  \class Profile
 *  \brief ABC for dielectric profiles
 *  \author Roberto Di Remigio
 *  \date 2014
 */

class Profile
{
public:
    Profile() {}
    virtual ~Profile() {}
    /*! The permittivity profile of the transition layer
     *  \param[out]  e the value of the dielectric constant at point r
     *  \param[out] de the value of the derivative of the dielectric constant
     *                 at point r
     *  \param[in]   r evaluation point
     */
    void operator()(double & e, double & de, const double r) const 
    {
        e = value(r);
        de = derivative(r);
    }
private:
    virtual double value(double point) const = 0;
    virtual double derivative(double point) const = 0;
};

#endif // PROFILE_HPP
