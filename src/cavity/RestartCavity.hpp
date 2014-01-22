#ifndef RESTARTCAVITY_HPP
#define RESTARTCAVITY_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include "Cavity.hpp"
#include "CavityData.hpp"
#include "CavityFactory.hpp"

/*! \file RestartCavity.hpp
 *  \class RestartCavity
 *  \brief A class for Restart cavity. 
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 */

class RestartCavity : public Cavity 
{
	public:
		RestartCavity(const std::string & _fname) : file(_fname) { makeCavity(); }
                virtual ~RestartCavity() {}
		virtual void makeCavity() { loadCavity(file); }
                friend std::ostream & operator<<(std::ostream & os, RestartCavity & cavity)
		{
			return cavity.printCavity(os);
		}
	private:
		std::string file;
                virtual std::ostream & printCavity(std::ostream & os);
};

namespace
{
	Cavity* createRestartCavity(const cavityData & _data)
	{
		return new RestartCavity(_data.filename);
	}
	const std::string RESTART("Restart");
	const bool registeredRestart = CavityFactory::TheCavityFactory().registerCavity(RESTART, createRestartCavity);
}

#endif // RESTARTCAVITY_HPP
