#ifndef SURFACEMAP_H
#define SURFACEMAP_H

#include <iostream>
#include <map>

#include "Config.h"

#include "SurfaceFunction.h"

class SurfaceMap
{
	public:
		typedef std::map<std::string, SurfaceFunction *> SurfaceFunctionMap;
		static SurfaceFunctionMap & TheMap()
		{
			static SurfaceFunctionMap obj;
			return obj;
		}

	private:
		SurfaceMap();
		SurfaceMap(const SurfaceMap & other);
		SurfaceMap & operator=(const SurfaceMap & other);
		~SurfaceMap();


};

#endif
