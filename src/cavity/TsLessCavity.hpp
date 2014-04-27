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

#ifndef TSLESSCAVITY_HPP
#define TSLESSCAVITY_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"


#include "Cavity.hpp"
#include "CavityData.hpp"
#include "CavityFactory.hpp"

/*! \file TsLessCavity.hpp
 *  \class TsLessCavity
 *  \brief A class for TsLess cavity.
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class is a wrapper for the Fortran routines generating
 *  a tessellationless grid on the cavity surface. The original
 *  algoritm is described in \cite Pomelli2004
 */

class TsLessCavity : public Cavity
{
public:
    TsLessCavity() {}
    TsLessCavity(const std::vector<Sphere> & _spheres, double _area,
                 double _probeRadius = 0.0,
                 double _minDistance = 0.1, int _derOrder = 4, double _minRadius = 100.0) :
        Cavity(_spheres), averageArea(_area), probeRadius(_probeRadius),
        minDistance(_minDistance),
        derOrder(_derOrder), minimalRadius(_minRadius) {
        build(10000, 200, 25000);
    }
    virtual ~TsLessCavity() {}
    friend std::ostream & operator<<(std::ostream & os, TsLessCavity & cavity) {
        return cavity.printCavity(os);
    }
private:
    double averageArea;
    double probeRadius;
    double minDistance;
    int derOrder;
    double minimalRadius;
    int addedSpheres;
    virtual std::ostream & printCavity(std::ostream & os);
    virtual void makeCavity() {
        build(10000, 200, 25000);
    }
    void build(int maxts, int maxsp, int maxvert);
};

namespace
{
    Cavity* createTsLessCavity(const cavityData & _data)
    {
        return new TsLessCavity(_data.spheres, _data.area, _data.probeRadius,
                                _data.minDistance, _data.derOrder, _data.minimalRadius);
    }
    const std::string TSLESS("TsLess");
    const bool registeredTsLess = CavityFactory::TheCavityFactory().registerCavity(
                                      TSLESS, createTsLessCavity);
}

#endif // TSLESSCAVITY_HPP
