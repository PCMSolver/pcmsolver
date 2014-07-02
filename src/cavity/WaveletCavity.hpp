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

#ifndef WAVELETCAVITY_HPP
#define WAVELETCAVITY_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include "EigenPimpl.hpp"

extern "C"
{
//#include "vector3.h"
}

#include "Vector3.hpp"

#include "Cavity.hpp"
#include "CavityData.hpp"
#include "CavityFactory.hpp"

/*! \file WaveletCavity.hpp
 *  \class WaveletCavity
 *  \brief A class for wavelet cavity.
 *  \author Luca Frediani
 *  \date 2011
 *
 *  This class is an interface to the C code wavcav for the generation
 *  of cavities according to the wavelet algorithms
 */

class WaveletCavity : public Cavity
{
public:
    WaveletCavity() {}
    WaveletCavity(const std::vector<Sphere> & _spheres, double _probeRadius,
                  int _patchLevel = 2, double _coarsity = 0.5) :
        Cavity(_spheres), probeRadius(_probeRadius), patchLevel(_patchLevel),
        coarsity(_coarsity) {
        uploadedDyadic = false;
        makeCavity();
    }
    virtual ~WaveletCavity() {};
    void readCavity(const std::string & filename);
    void uploadPoints(int quadLevel, vector3 **** T_, bool isPWL);
    unsigned int getNPatches() { return nPatches; }
    unsigned int getNPatches() const { return nPatches; }
    unsigned int getNLevels() { return nLevels; }
    unsigned int getNLevels() const { return nLevels; }
    unsigned int getNPoints() { return nPoints; }
    unsigned int getNPoints() const { return nPoints; }
    Eigen::Vector3d getNodePoint(int i) { return nodePoint[i]; }
    Eigen::Vector3d getNodePoint(int i) const { return nodePoint[i]; }
    Eigen::Vector3i getNodeIndex(int i) { return nodeIndex[i]; }
    Eigen::Vector3i getNodeIndex(int i) const { return nodeIndex[i]; }
    friend std::ostream & operator<<(std::ostream & os, WaveletCavity & cavity) {
        return cavity.printCavity(os);
    }
    void compFakePotential();
private:
    double probeRadius;
    int patchLevel;
    double coarsity;
    void uploadPointsPWC(int quadLevel, vector3 **** T_);
    void uploadPointsPWL(int quadLevel, vector3 **** T_);
    std::vector<Eigen::Vector3d> nodePoint;
    std::vector<Eigen::Vector3i> nodeIndex;
    unsigned int nPatches;
    unsigned int nLevels;
    unsigned int nPoints;
    bool uploadedDyadic;
    void writeInput(std::string &fileName);
    virtual std::ostream & printCavity(std::ostream & os);
    virtual void makeCavity();
};

namespace
{
    Cavity* createWaveletCavity(const cavityData & _data)
    {
        return new WaveletCavity(_data.spheres, _data.probeRadius, _data.patchLevel,
                                 _data.coarsity);
    }
    const std::string WAVELET("Wavelet");
    const bool registeredWavelet = CavityFactory::TheCavityFactory().registerCavity(
                                       WAVELET, createWaveletCavity);
}

#endif // WAVELETCAVITY_HPP
