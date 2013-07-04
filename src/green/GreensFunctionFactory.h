#ifndef GREENSFUNCTIONFACTORY_H
#define GREENSFUNCTIONFACTORY_H

#include <iostream>
#include <string>
#include <map>

#include "GreensFunction.h"

/*!
 *	\file GreensFunctionFactory.h
 *	\class GreensFunctionFactory
 *	\brief Implementation of the Factory Method for Green's functions. 
 *	\author Roberto Di Remigio
 *	\date 2013 
 *
 * 	Factory method implementation shamelessly copied from "Modern C++ Design" of A. Alexandrescu.
 * 	It is implemented as a Singleton.
 */

class GreensFunctionFactory
{
	public:
		/*!
		 * Callback function for Green's function creation.
		 */
		typedef GreensFunction * (*createGreensFunctionCallback)(int how_, double epsilon_, double kappa_);
	private:
		/*!
		 * A map from the Green's function type identifier (a string) to its callback function.
		 */
		typedef std::map<std::string, createGreensFunctionCallback> CallbackMap;
	public:
		/*!
		 * \brief Returns true if registration of the greenID was successful
		 * \param greenID the Green's function identification string
		 * \param createFunction the creation function related to the Green's function type given
		 */
		bool registerGreensFunction(const std::string & greenID, createGreensFunctionCallback createFunction);
		/*!
		 * \brief Returns true if greenID was already registered
		 * \param greenID the Green's function identification string
		 */
		bool unRegisterGreensFunction(const std::string & greenID);
		/*! 
		 * Calls the appropriate creation function, based on the passed greenID
		 */
		GreensFunction * createGreensFunction(const std::string & greenID, int how_, double epsilon_ = 1.0, double kappa_ = 0.0);
		/*!
		 * Unique point of access to the unique instance of the GreensFunctionFactory
		 */
		static GreensFunctionFactory& TheGreensFunctionFactory() 
		{
			static GreensFunctionFactory obj;
			return obj;
		}
	private:
		GreensFunctionFactory(){}
		/// Copy constructor is made private
		GreensFunctionFactory(const GreensFunctionFactory &other);
		GreensFunctionFactory& operator=(const GreensFunctionFactory &other);
        	~GreensFunctionFactory(){}
		CallbackMap callbacks;	
};

#endif // GREENSFUNCTIONFACTORY_H
