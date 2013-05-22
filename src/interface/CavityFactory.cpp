#include <stdexcept>

#include "CavityFactory.h"

;

bool CavityFactory::RegisterCavity(int cavityID, CreateCavityCallback createFunction)
{
	return callbacks.insert(CallbackMap::value_type(cavityID, createFunction)).second;
}

bool CavityFactory::UnregisterCavity(int cavityID)
{
	return callbacks.erase(cavityID) == 1;
}

Cavity * CavityFactory::CreateCavity(int cavityID)
{
	CallbackMap::const_iterator i = callbacks.find(cavityID);
	if (i == callbacks.end()) 
	{
		// The cavityID was not found
                throw std::runtime_error("Unknown cavity ID.");
	}
	// Invoke the creation function
	return (i->second)();
}
