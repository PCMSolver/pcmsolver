#ifndef CAVITYFACTORY_H
#define CAVITYFACTORY_H

#include <iostream>
#include <string>
#include <map>

#include "Config.h"

#include <Eigen/Dense>

#include "Cavity.h"
#include "CavityOfSpheres.h"
#include "GePolCavity.h"
#include "WaveletCavity.h"

/*
 * Factory method implementation shamelessly copied from
 * "Modern C++ Design" of A. Alexandrescu
 */


class CavityFactory {
	public:
		typedef Cavity * (*CreateCavityCallback)();
	private:
		typedef std::map<int, CreateCavityCallback> CallbackMap;
	public:
		// Returns 'true' if registration was successful
		bool RegisterCavity(int cavityID, CreateCavityCallback createFunction);
		// Returns 'true' if registration it the cavityType was registered before
		bool UnregisterCavity(int cavityID);
		Cavity * CreateCavity(int cavityID);
		static CavityFactory& CreateCavityFactory() 
		{
			static CavityFactory obj;
			return obj;
		}
	private:
		CavityFactory();
		/// Copy constructor is made private
		CavityFactory(const CavityFactory &other);
		CavityFactory& operator=(const CavityFactory &other);
        	~CavityFactory(){}
		CallbackMap callbacks;	
}

#endif
