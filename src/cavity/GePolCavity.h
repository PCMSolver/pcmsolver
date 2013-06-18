#ifndef GEPOLCAVITY
#define GEPOLCAVITY

#include <iostream>
#include <string>
#include <vector>

#include <Config.h>

#include "Cavity.h"
#include "CavityFactory.h"

/*

  C++ inteface and wrapper for GePol cavity class
  written by Krzysztof Mozgawa, 2011

*/

class GePolCavity : public Cavity 
{
 public:
    GePolCavity(){}
    GePolCavity(const std::vector<Sphere> & _spheres, double _area, double _probeRadius = 0.0, bool _addSpheres = false) : 
	Cavity(_spheres), averageArea(_area), probeRadius(_probeRadius), addSpheres(_addSpheres)  
           {
		   setMaxAddedSpheres(addSpheres);
		   makeCavity(10000, 10000000);
           }
    ~GePolCavity(){}
    void makeCavity(int maxts, int lwork);
    void makeCavity();
    void writeOutput(std::string &filename);
    void setMaxAddedSpheres(bool add = true, int maxAdd = 100);
    double getProbeRadius() { return probeRadius; }
    void setProbeRadius( double probeRadius );

 protected:
    virtual std::ostream & printObject(std::ostream & os);
        
 private:
    double averageArea;
    bool addSpheres;
    double probeRadius;
    int maxAddedSpheres;
    int addedSpheres;
};

namespace
{
	Cavity* createGePolCavity(const std::vector<Sphere> & _spheres, double _area, double _probeRadius = 0.0, 
			bool _addSpheres = false, int _patchLevel = 2, double _coarsity = 0.5)
	{
		return new GePolCavity(_spheres, _area, _probeRadius, _addSpheres);
	}
	const std::string GEPOL("GePol");
	const bool registeredGePol = CavityFactory::TheCavityFactory().registerCavity(GEPOL, createGePolCavity);
}

#endif
