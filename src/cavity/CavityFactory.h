#ifndef CAVITYFACTORY_H
#define CAVITYFACTORY_H

#include <iostream>
#include <string>
#include <map>

#include "Config.h"

#include <Eigen/Dense>

#include "Cavity.h"

/*
 * Factory method implementation shamelessly copied from
 * "Modern C++ Design" of A. Alexandrescu
 */


class CavityFactory {
	public:
		typedef Cavity * (*CreateCavityCallback)(const std::vector<Sphere> & _spheres, double _area, double _probeRadius, 
							 bool _addSpheres, int _patchLevel, double _coarsity);
	private:
		typedef std::map<std::string, CreateCavityCallback> CallbackMap;
	public:
		// Returns 'true' if registration was successful
		bool registerCavity(std::string cavityID, CreateCavityCallback createFunction);
		// Returns 'true' if registration it the cavityType was registered before
		bool unRegisterCavity(std::string cavityID);
		Cavity * CreateCavity(std::string cavityID, const std::vector<Sphere> & _spheres, double _area, double _probeRadius = 0.0, 
				      bool _addSpheres = false, int _patchLevel = 2, double _coarsity = 0.5);
		static CavityFactory& TheCavityFactory() 
		{
			static CavityFactory obj;
			return obj;
		}
	private:
		CavityFactory(){}
		/// Copy constructor is made private
		CavityFactory(const CavityFactory &other);
		CavityFactory& operator=(const CavityFactory &other);
        	~CavityFactory(){}
		CallbackMap callbacks;	
};

#endif
