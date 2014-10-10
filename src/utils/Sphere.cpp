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

#include "Sphere.hpp"

#include <ostream>
#include <stdexcept>

#include "Config.hpp"


std::ostream & Sphere::printObject(std::ostream & os)
{
    os << "Sphere radius " << radius_ << std::endl;
    os << "Sphere center\n" << center_;

    return os;
}

inline void swap(Sphere & left, Sphere & right)
{
    using std::swap;
    swap(left.center_, right.center_);
    swap(left.radius_, right.radius_);
    swap(left.colour_, right.colour_);
}

inline void Sphere::swap(Sphere & other)
{
    using std::swap;
    swap(this->center_, other.center_);
    swap(this->radius_, other.radius_);
    swap(this->colour_, other.colour_);
}

Sphere & Sphere::operator=(Sphere other)
{
    if (this == &other) // Check for self-assignment
        throw std::runtime_error("Are you trying to self-assign?");

    ::swap(*this, other);
    return *this;
}
