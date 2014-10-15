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

#ifndef WAVELETCAVITY_HPP
#define WAVELETCAVITY_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

extern "C"
{
//#include "vector3.h"
}

#include "Vector3.hpp"
#include "Interpolation.hpp"

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
    WaveletCavity(const std::vector<Sphere> & s, double pr, int pl = 2, double c = 0.5) :
        Cavity(s), probeRadius_(pr), patchLevel_(pl), coarsity_(c) {
        uploadedDyadic_ = false;
        makeCavity();
    }
    virtual ~WaveletCavity() {};
    void readCavity(const std::string & filename);
    void uploadPoints(int quadLevel, Interpolation *interp);
    unsigned int getNPatches() { return nPatches_; }
    unsigned int getNPatches() const { return nPatches_; }
    unsigned int getNLevels() { return nLevels_; }
    unsigned int getNLevels() const { return nLevels_; }
    unsigned int getNPoints() { return nPoints_; }
    unsigned int getNPoints() const { return nPoints_; }
    Eigen::Vector3d getNodePoint(int i) { return nodePoint_[i]; }
    Eigen::Vector3d getNodePoint(int i) const { return nodePoint_[i]; }
    Eigen::Vector3i getNodeIndex(int i) { return nodeIndex_[i]; }
    Eigen::Vector3i getNodeIndex(int i) const { return nodeIndex_[i]; }
    friend std::ostream & operator<<(std::ostream & os, WaveletCavity & cavity) {
        return cavity.printCavity(os);
    }
    void compFakePotential();
private:
    double probeRadius_;
    int patchLevel_;
    double coarsity_;
    void uploadPointsPWC(int quadLevel, Interpolation *interp);
    void uploadPointsPWL(int quadLevel, Interpolation *interp);
    std::vector<Eigen::Vector3d> nodePoint_;
    std::vector<Eigen::Vector3i> nodeIndex_;
    unsigned int nPatches_;
    unsigned int nLevels_;
    unsigned int nPoints_;
    bool uploadedDyadic_;
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
    const std::string WAVELET("WAVELET");
    const bool registeredWavelet = CavityFactory::TheCavityFactory().registerCavity(
                                       WAVELET, createWaveletCavity);
}

#endif // WAVELETCAVITY_HPP
