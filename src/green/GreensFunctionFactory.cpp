#include <stdexcept>

#include "GreensFunctionFactory.h"

template <typename T>
bool GreensFunctionFactory<T>::registerGreensFunction(std::string greenID, createGreensFunctionCallback createFunction)
{
	return callbacks.insert(CallbackMap::value_type(greenID, createFunction)).second;
}

template <typename T>
bool GreensFunctionFactory<T>::unRegisterGreensFunction(std::string greenID)
{
	return callbacks.erase(greenID) == 1;
}

template <typename T>
GreensFunction<T> * GreensFunctionFactory<T>::createGreensFunction(std::string greenID) 
{
	CallbackMap::const_iterator i = callbacks.find(greenID);
	if (i == callbacks.end()) 
	{
		// The cavityID was not found
                throw std::runtime_error("Unknown Green's Function ID.");
	}
	// Invoke the creation function
	return (i->second)();
}
