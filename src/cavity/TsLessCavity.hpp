/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2016 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */
/* pcmsolver_copyright_end */

#ifndef TSLESSCAVITY_HPP
#define TSLESSCAVITY_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"

#include "Cavity.hpp"

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
    TsLessCavity(const Molecule & molec, double a, double pr, double minR, double minD, int der) :
        Cavity(molec), averageArea_(a), probeRadius_(pr), minimalRadius_(minR), minDistance_(minD), derOrder_(der)
    {
	    std::string checkpointName = "TsLessCavity::build";
	    TIMER_ON(checkpointName);
        build(10000, 200, 25000);
        TIMER_OFF(checkpointName);
    }
    TsLessCavity(const std::vector<Sphere> & sph, double a, double pr, double minR, double minD, int der) :
        Cavity(sph), averageArea_(a), probeRadius_(pr), minimalRadius_(minR), minDistance_(minD), derOrder_(der)
	{
	  std::string checkpointName = "TsLessCavity::build";
	  TIMER_ON(checkpointName);
	  build(10000, 200, 25000);
	  TIMER_OFF(checkpointName);
	}
    virtual ~TsLessCavity() {}
    friend std::ostream & operator<<(std::ostream & os, TsLessCavity & cavity) {
        return cavity.printCavity(os);
    }
private:
    double averageArea_;
    double probeRadius_;
    double minimalRadius_;
    double minDistance_;
    int derOrder_;
    int addedSpheres;
    virtual std::ostream & printCavity(std::ostream & os);
    virtual void makeCavity() { build(10000, 200, 25000); }
    /*! \brief Driver for TsLess Fortran module.
     *  \param[in]   maxts maximum number of tesserae
     *  \param[in]   maxsp maximum number of spheres (original + added)
     *  \param[in] maxvert maximum number of vertices
     */
    void build(size_t maxts, size_t maxsp, size_t maxvert);
};

#endif // TSLESSCAVITY_HPP
