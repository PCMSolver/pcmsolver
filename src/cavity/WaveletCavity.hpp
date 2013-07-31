#ifndef WAVELETCAVITY_HPP
#define WAVELETCAVITY_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

extern "C"
{
#include "vector3.h"
}

#include "Cavity.hpp"
#include "CavityFactory.hpp"

class WaveletCavity : public Cavity 
{
	public:
		WaveletCavity(){}
                WaveletCavity(const std::vector<Sphere> & _spheres, double _probeRadius, int _patchLevel = 2, double _coarsity = 0.5) :
            	   Cavity(_spheres), probeRadius(_probeRadius), patchLevel(_patchLevel), coarsity(_coarsity) 
                       {
			       uploadedDyadic = false;
            		       makeCavity();
                       }
                virtual ~WaveletCavity(){};
                void makeCavity();
                void readCavity(const std::string & filename);
                void uploadPoints(int quadLevel, vector3 **** T_, bool isPWL);
                unsigned int getNPatches() { return nPatches; }
                unsigned int getNLevels() { return nLevels; }
                unsigned int getNPoints() { return nPoints; }
                Eigen::Vector3d getNodePoint(int i) { return nodePoint[i]; }
                Eigen::Vector3i getNodeIndex(int i) { return nodeIndex[i]; }
                friend std::ostream & operator<<(std::ostream & os, WaveletCavity & cavity)
		{
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
};

namespace
{
	Cavity* createWaveletCavity(const std::vector<Sphere> & _spheres, double _area, double _probeRadius = 0.0, 
		    bool _addSpheres = false, int _patchLevel = 2, double _coarsity = 0.5)
	{
		return new WaveletCavity(_spheres, _probeRadius, _patchLevel, _coarsity);
        }
	const std::string WAVELET("Wavelet");
	const bool registeredWavelet = CavityFactory::TheCavityFactory().registerCavity(WAVELET, createWaveletCavity);
}

#endif // WAVELETCAVITY_HPP
