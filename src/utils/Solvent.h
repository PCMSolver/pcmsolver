#ifndef SOLVENT_H
#define SOLVENT_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <Config.h>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/*!
 * \brief Class describing a solvent
 * \author Roberto Di Remigio
 * \date 2011
 *
 * A Solvent object contains all the solvent-related experimental data
 * needed to set up the Green's functions and the non-electrostatic
 * terms calculations.
 * 
 * These data are taken from the DALTON2011 internal implementation of
 * the Polarizable Continuum Model.
 */

class Solvent;

typedef std::map< std::string, Solvent > SolventMap;

class Solvent {
 public:
    Solvent(){}
    Solvent( const string & _name, double _epsStatic, 
             double _epsOptical, double _radius );
    ~Solvent(){}
    const string getName(){ return name; }
    double getEpsStatic(){ return epsStatic; }
    double getEpsOptical(){ return epsOptical; }
    double getRadius() const { return probeRadius; }
    void setRadius(double r){ probeRadius = r; }
    static SolventMap initSolventMap();
    friend std::ostream & operator<<(std::ostream & os, Solvent & obj) {
        return obj.printObject(os);
    }
 private:
    ostream & printObject(ostream & os);
    string name;
    double epsStatic;
    double epsOptical;
    double probeRadius;
};

#endif
