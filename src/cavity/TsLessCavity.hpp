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
 */

class TsLessCavity : public Cavity 
{
	public:
		TsLessCavity(){}
                TsLessCavity(const std::vector<Sphere> & _spheres, double _area, double _probeRadius = 0.0, 
			     double _minDistance = 0.1, int _derOrder = 4, double _minRadius = 100.0) :  
                    Cavity(_spheres), averageArea(_area), probeRadius(_probeRadius), minDistance(_minDistance), 
		    derOrder(_derOrder), minimalRadius(_minRadius)
                       {
                    	   makeCavity(10000, 10000000);
                       }
                virtual ~TsLessCavity() {}
                void makeCavity(int maxts, int lwork);
                void makeCavity();
                friend std::ostream & operator<<(std::ostream & os, TsLessCavity & cavity)
		{
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
};

namespace
{
	Cavity* createTsLessCavity(const cavityData & _data)
	{
		return new TsLessCavity(_data.spheres, _data.area, _data.probeRadius, _data.minDistance, _data.derOrder, _data.minimalRadius);
	}
	const std::string TSLESS("TsLess");
	const bool registeredTsLess = CavityFactory::TheCavityFactory().registerCavity(TSLESS, createTsLessCavity);
}

#endif // TSLESSCAVITY_HPP
