#ifndef SOLVENT_H
#define SOLVENT_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/*

  A Solvent class
  written by Roberto Di Remigio, 2011

  "Why use a structure when you can define a class?"
  
  vector<Solvent> SolventData[] should contain all the solvent-related
  data needed to set up the Green's functions and the non-electrostatic
  terms calculations.
  
  These data are taken from the DALTON2011 internal implementation of
  the Polarizable Continuum Model.

*/

class Solvent;

typedef std::map< std::string, Solvent > SolventMap;

class Solvent {
 public:
    Solvent(){}
    Solvent( const string & name, double epsstatic, 
             double epsoptical, double radius );
    ~Solvent(){}
    const string getName(){ return name; }
    double getEpsStatic(){ return epsStatic; }
    double getEpsOptical(){ return epsOptical; }
    const double getRadius() const { return probeRadius; }
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
