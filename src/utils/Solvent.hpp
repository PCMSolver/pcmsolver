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
    typedef std::map<std::string, Solvent> SolventMap;

    Solvent() {}
    Solvent(const std::string & name, double epsStatic, double epsOptical,
            double radius )
        : name_(name), epsStatic_(epsStatic), epsOptical_(epsOptical),
          probeRadius_(radius) {}
    ~Solvent() {}

    std::string name() const { return name_; }
    double epsStatic() const { return epsStatic_; }
    double epsOptical() const { return epsOptical_; }
    double probeRadius() const { return (probeRadius_ / convertBohrToAngstrom); }
    void probeRadius(double radius) { probeRadius_ = radius; }

    /*! \brief Returns the map between solvent names and Solvent objects.
     *
     *  This map contains solvent data taken from the DALTON2011 internal
     *  implementation of the Polarizable Continuum Model.
     */
    static SolventMap & initSolventMap();
    friend std::ostream & operator<<(std::ostream & os, Solvent & solvent) {
        return solvent.printSolvent(os);
    }
private:
    std::string name_;
    double epsStatic_;
    double epsOptical_;
    double probeRadius_;
    std::ostream & printSolvent(std::ostream & os);
};

#endif // SOLVENT_HPP
