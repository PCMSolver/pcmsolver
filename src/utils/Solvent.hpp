#ifndef SOLVENT_HPP
#define SOLVENT_HPP

#include <iosfwd>
#include <string>
#include <map>

#include "Config.hpp"

#include "PhysicalConstants.hpp"

/*! \file Solvent.hpp
 *  \class Solvent
 *  \brief Class describing a solvent.
 *  \author Roberto Di Remigio
 *  \date 2011
 *
 * A Solvent object contains all the solvent-related experimental data
 * needed to set up the Green's functions and the non-electrostatic
 * terms calculations. 
 */

class Solvent 
{
	public:
		/*! \brief typedef for the map between solvent name and Solvent object.
		 */
		typedef std::map< std::string, Solvent > SolventMap;

		Solvent(){}
                Solvent(const std::string & _name, double _epsStatic, double _epsOptical, double _radius )
			: name(_name), epsStatic(_epsStatic), epsOptical(_epsOptical), probeRadius(_radius) {}
                ~Solvent(){}

                std::string getName() const { return name; }
                double getEpsStatic() const { return epsStatic; }
                double getEpsOptical() const { return epsOptical; }
                double getRadius() const { return (probeRadius / convertBohrToAngstrom); }
                void setRadius(double _radius) { probeRadius = _radius; }

		/*! \brief Returns the map between solvent names and Solvent objects.
		 *
		 *  This map contains solvent data taken from the DALTON2011 internal 
		 *  implementation of the Polarizable Continuum Model.
		 */
                static SolventMap & initSolventMap();
		friend std::ostream & operator<<(std::ostream & os, Solvent & solvent) 
		{
                    return solvent.printSolvent(os);
                }
	private:
		std::string name;        
                double epsStatic;
                double epsOptical;
                double probeRadius;
		std::ostream & printSolvent(std::ostream & os);
};

#endif // SOLVENT_HPP
