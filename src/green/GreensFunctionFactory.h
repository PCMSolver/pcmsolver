#ifndef GREENSFUNCTIONFACTORY_H
#define GREENSFUNCTIONFACTORY_H

#include <iostream>
#include <map>

//#include "GreensFunctionInterface.h"
#include "GreensFunction.h"

/*! \file GreensFunctionFactory.h
 *  \class GreensFunctionFactory
 *  \brief A templatized factory for Green's Functions.
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This factory is a Singleton.
 */

template <typename T>
class GreensFunctionFactory
{
	public:
		/*!
		 * Callback function for Green's function creation.
		 */
		typedef GreensFunction<T> * (*createGreensFunctionCallback)();
	private:
		/*!
		 * A map from the Green's function derivative strategy identifier (a string) to its callback function.
		 */
		typedef std::map<std::string, createGreensFunctionCallback> CallbackMap;
	public:
		/*!
		 * \brief Returns true if registration of the greenDerID was successful
		 * \param greenID the Green's function identification string
		 * \param createFunction the creation function related to the identification string given
		 */
		bool registerGreensFunction(std::string greenID, createGreensFunctionCallback createFunction);
		/*!
		 * \brief Returns true if greenDerID was already registered
		 * \param[in] greenID the Green's function identification string
		 */
		bool unRegisterGreensFunction(std::string greenID);
		/*! 
		 * Calls the appropriate creation function, based on the passed greenID
		 */
		GreensFunction<T> * createGreensFunction(std::string greenDerID);
		/*!
		 * Unique point of access to the unique instance of the GreensFunctionFactory<T>
		 */
		static GreensFunctionFactory<T>& TheGreensFunctionFactory() 
		{
			static GreensFunctionFactory obj;
			return obj;
		}
	
	private:
		GreensFunctionFactory<T>(){}
		~GreensFunctionFactory<T>(){}
		GreensFunctionFactory<T>(const GreensFunctionFactory<T> &other);
		GreensFunctionFactory<T> & operator=(const GreensFunctionFactory<T> &other);
		CallbackMap callbacks;
};

#endif // GREENSFUNCTIONFACTORY_H
