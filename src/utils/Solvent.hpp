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

#ifndef SOLVENT_HPP
#define SOLVENT_HPP

#include <iosfwd>
#include <string>
#include <map>

#include "Config.hpp"

/*! \file Solvent.hpp
 *  \class Solvent
 *  \brief Class describing a solvent.
 *  \author Roberto Di Remigio
 *  \date 2011
 *
 * A Solvent object contains all the solvent-related experimental data
 * needed to set up the Green's functions and the non-electrostatic
 * terms calculations.
 */

class Solvent
{
public:
    /*! \brief typedef for the map between solvent name and Solvent object.
     */
    typedef std::map<std::string, Solvent> SolventMap;

    Solvent() {}
    Solvent(const std::string & name, double epsStatic, double epsOptical,
            double radius )
        : name_(name), epsStatic_(epsStatic), epsOptical_(epsOptical),
          probeRadius_(radius) {}
    ~Solvent() {}

    std::string name() const { return name_; }
    double epsStatic() const { return epsStatic_; }
    double epsOptical() const { return epsOptical_; }
    // The built-in probe radii are in Angstrom
    double probeRadius() const { return probeRadius_; }
    void probeRadius(double radius) { probeRadius_ = radius; }

    /*! \brief Returns the map between solvent names and Solvent objects.
     *
     *  This map contains solvent data taken from the DALTON2011 internal
     *  implementation of the Polarizable Continuum Model.
     */
    static SolventMap & initSolventMap();
    friend std::ostream & operator<<(std::ostream & os, Solvent & solvent) {
        return solvent.printSolvent(os);
    }
private:
    std::string name_;
    double epsStatic_;
    double epsOptical_;
    double probeRadius_;
    std::ostream & printSolvent(std::ostream & os);
};

#endif // SOLVENT_HPP
