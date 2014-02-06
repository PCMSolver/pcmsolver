#ifndef GEPOLCAVITY_HPP
#define GEPOLCAVITY_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"

#include "Cavity.hpp"
#include "CavityData.hpp"
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
	private:
		enum pGroup { C1, Cs, C2, Ci, C2h, D2, C2v, D2h };
	public:
		GePolCavity(){}
                GePolCavity(const std::vector<Sphere> & _spheres, double _area, double _probeRadius = 0.0, double _minRadius = 100.0, int _pGroup = C1) :  
                    Cavity(_spheres), averageArea(_area), probeRadius(_probeRadius), minimalRadius(_minRadius), pointGroup(_pGroup) 
                       {
                    	   build(10000, 200, 25000);
                       }
                virtual ~GePolCavity(){}
                friend std::ostream & operator<<(std::ostream & os, GePolCavity & cavity)
		{
			return cavity.printCavity(os);
		}
	private:
                double averageArea;  
                double probeRadius;
		double minimalRadius;
		int pointGroup;
                int addedSpheres;
                virtual std::ostream & printCavity(std::ostream & os);
                virtual void makeCavity()
		{
		        build(10000, 200, 25000);
		}
	        /*! \brief Driver for PEDRA Fortran module. 
		 *  \param[in] maxts maximum number of tesserae
		 *  \param[in] maxsp maximum number of spheres (original + added)
		 *  \param[in] maxvert maximum number of vertices
		 */
                void build(int maxts, int maxsp, int maxvert);
};

namespace
{
	Cavity* createGePolCavity(const cavityData & _data)
	{
		return new GePolCavity(_data.spheres, _data.area, _data.probeRadius, _data.minimalRadius);
	}
	const std::string GEPOL("GePol");
	const bool registeredGePol = CavityFactory::TheCavityFactory().registerCavity(GEPOL, createGePolCavity);
}

#endif // GEPOLCAVITY_HPP
