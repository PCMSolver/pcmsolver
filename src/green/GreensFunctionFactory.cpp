#include <stdexcept>

#include "GreensFunctionFactory.h"

bool GreensFunctionFactory::registerGreensFunction(const std::string & greenID, createGreensFunctionCallback createFunction)
{
	return callbacks.insert(CallbackMap::value_type(greenID, createFunction)).second;
}

bool GreensFunctionFactory::unRegisterGreensFunction(const std::string & greenID)
{
	return callbacks.erase(greenID) == 1;
}

GreensFunction * GreensFunctionFactory::createGreensFunction(const std::string & greenID, int how_, double epsilon_, double kappa_) 
{
	CallbackMap::const_iterator i = callbacks.find(greenID);
	if (i == callbacks.end()) 
	{
		// The greenID was not found
                throw std::runtime_error("The unknown Green's function ID " + greenID + " occurred in GreensFunctionFactory.");
	}
	// Invoke the creation function
	return (i->second)(how_, epsilon_, kappa_);
}
