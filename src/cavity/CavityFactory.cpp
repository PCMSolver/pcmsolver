#include "CavityFactory.hpp"

#include <string>
#include <stdexcept>
#include <vector>

#include "Config.hpp"

#include "Sphere.hpp"

bool CavityFactory::registerCavity(std::string cavityID, createCavityCallback createFunction)
{
	return callbacks.insert(CallbackMap::value_type(cavityID, createFunction)).second;
}

bool CavityFactory::unRegisterCavity(std::string cavityID)
{
	return callbacks.erase(cavityID) == 1;
}

Cavity * CavityFactory::createCavity(std::string cavityID, const std::vector<Sphere> & _spheres, double _area, double _probeRadius, 
				     bool _addSpheres, int _patchLevel, double _coarsity)
{
	CallbackMap::const_iterator i = callbacks.find(cavityID);
	if (i == callbacks.end()) 
	{
		// The cavityID was not found
                throw std::runtime_error("The unknown cavity ID " + cavityID + " occurred in CavityFactory.");
	}
	// Invoke the creation function
	return (i->second)(_spheres, _area, _probeRadius, _addSpheres, _patchLevel, _coarsity);
}
