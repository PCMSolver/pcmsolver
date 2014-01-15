#ifndef GEPOLCAVITY_HPP
#define GEPOLCAVITY_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"

#include "Cavity.hpp"
#include "CavityFactory.hpp"

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
		GePolCavity(){}
                GePolCavity(const std::vector<Sphere> & _spheres, double _area, double _probeRadius = 0.0, bool _addSpheres = false) :  
                    Cavity(_spheres), averageArea(_area), probeRadius(_probeRadius), addSpheres(_addSpheres), maxAddedSpheres(100) 
                       {
                    	   makeCavity(10000, 10000000);
                       }
                virtual ~GePolCavity(){}
                void makeCavity(int maxts, int lwork);
                void makeCavity();
                double getProbeRadius() { return probeRadius; }
                void setProbeRadius( double probeRadius );
                friend std::ostream & operator<<(std::ostream & os, GePolCavity & cavity)
		{
			return cavity.printCavity(os);
		}
	private:
                double averageArea;  
                double probeRadius;
                bool addSpheres;
                int maxAddedSpheres;
                int addedSpheres;
                virtual std::ostream & printCavity(std::ostream & os);  
};

namespace
{
	Cavity* createGePolCavity(const cavityData & _data)
	{
		return new GePolCavity(_data.spheres, _data.area, _data.probeRadius, _data.addSpheres);
	}
	const std::string GEPOL("GePol");
	const bool registeredGePol = CavityFactory::TheCavityFactory().registerCavity(GEPOL, createGePolCavity);
}

#endif // GEPOLCAVITY_HPP
