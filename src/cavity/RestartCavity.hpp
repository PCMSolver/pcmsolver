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

#ifndef RESTARTCAVITY_HPP
#define RESTARTCAVITY_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include "Cavity.hpp"

/*! \file RestartCavity.hpp
 *  \class RestartCavity
 *  \brief A class for Restart cavity.
 *  \author Roberto Di Remigio
 *  \date 2014
 */

class RestartCavity __final : public Cavity
{
public:
    RestartCavity(const std::string & _fname) : file(_fname) {
        makeCavity();
    }
    virtual ~RestartCavity() {}
    virtual void makeCavity() __override {
        loadCavity(file);
    }
    friend std::ostream & operator<<(std::ostream & os, RestartCavity & cavity) {
        return cavity.printCavity(os);
    }
private:
    std::string file;
    virtual std::ostream & printCavity(std::ostream & os) __override;
};

#endif // RESTARTCAVITY_HPP
