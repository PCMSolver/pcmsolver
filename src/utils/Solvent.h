#ifndef SOLVENT_H
#define SOLVENT_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <Config.h>

/*!
 * \brief Class describing a solvent
 * \author Roberto Di Remigio
 * \date 2011
 *
 * A Solvent object contains all the solvent-related experimental data
 * needed to set up the Green's functions and the non-electrostatic
 * terms calculations. 
 */

class Solvent 
{
	public:
		/*!
		 * typedef for the map between solvent name and Solvent object.
		 */
		typedef std::map< std::string, Solvent > SolventMap;

		Solvent(){}
                Solvent(const std::string & _name, double _epsStatic, double _epsOptical, double _radius )
			: name(_name), epsStatic(_epsStatic), epsOptical(_epsOptical), probeRadius(_radius) {}
                ~Solvent(){}

                const std::string getName(){ return name; }
                double getEpsStatic(){ return epsStatic; }
                double getEpsOptical(){ return epsOptical; }
                double getRadius() const { return probeRadius; }
                void setRadius(double _radius){ probeRadius = _radius; }

		/*!
		 * This map contains solvent data taken from the DALTON2011 internal 
		 * implementation of the Polarizable Continuum Model.
		 */
                static SolventMap & initSolventMap();
                
		friend std::ostream & operator<<(std::ostream & os, Solvent & obj) 
		{
                    return obj.printObject(os);
                }

	private:
		std::string name;        
                double epsStatic;
                double epsOptical;
                double probeRadius;
		std::ostream & printObject(std::ostream & os);
};

#endif
