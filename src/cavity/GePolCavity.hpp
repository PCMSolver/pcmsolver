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

#ifndef GEPOLCAVITY_HPP
#define GEPOLCAVITY_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"

#include <boost/current_function.hpp>

#include "Cavity.hpp"
#include "CavityData.hpp"
#include "CavityFactory.hpp"
#include "TimerInterface.hpp"

/*! \file GePolCavity.hpp
 *  \class GePolCavity
 *  \brief A class for GePol cavity.
 *  \author Krzysztof Mozgawa
 *  \date 2011
 *
 *  This class is an interface to the Fortran code PEDRA for the generation
 *  of cavities according to the GePol algorithm.
 */

class GePolCavity : public Cavity
{
public:
    GePolCavity() {}
    GePolCavity(const Molecule & molec, double a, double pr, double minR) :
        Cavity(molec), averageArea(a), probeRadius(pr), minimalRadius(minR)
	{
	  std::string checkpointName = "GePolCavity::build";
	  TIMER_ON(checkpointName);
	  build(10000, 200, 25000);
	  TIMER_OFF(checkpointName);
	}
    GePolCavity(const std::vector<Sphere> & sph, double a, double pr, double minR) :
	Cavity(sph), averageArea(a), probeRadius(pr), minimalRadius(minR)
	{
	  std::string checkpointName = "GePolCavity::build";
	  TIMER_ON(checkpointName);
	  build(10000, 200, 25000);
	  TIMER_OFF(checkpointName);
	}
    virtual ~GePolCavity() {}
    friend std::ostream & operator<<(std::ostream & os, GePolCavity & cavity) {
        return cavity.printCavity(os);
    }
private:
    double averageArea;
    double probeRadius;
    double minimalRadius;
    int addedSpheres;
    virtual std::ostream & printCavity(std::ostream & os);
    virtual void makeCavity() { build(10000, 200, 25000); }
    /*! \brief Driver for PEDRA Fortran module.
     *  \param[in]   maxts maximum number of tesserae
     *  \param[in]   maxsp maximum number of spheres (original + added)
     *  \param[in] maxvert maximum number of vertices
     */
    void build(int maxts, int maxsp, int maxvert);
};

namespace
{
    Cavity* createGePolCavity(const cavityData & data)
    {
        return new GePolCavity(data.molecule, data.area, data.probeRadius,
                               data.minimalRadius);
    }
    const std::string GEPOL("GEPOL");
    const bool registeredGePol = CavityFactory::TheCavityFactory().registerCavity(GEPOL,
                                 createGePolCavity);
}

#endif // GEPOLCAVITY_HPP
