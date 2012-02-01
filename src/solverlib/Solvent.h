#ifndef SOLVENT_H
#define SOLVENT_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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

class Solvent {
 public:
    Solvent(){}
    Solvent( const string & name, double epsstatic, 
             double epsoptical, double radius );
    ~Solvent(){}
    string getSolventName(){ return solventName; }
    void setSolventName( const string & name ){ solventName = name; }
    double getSolventEpsStatic(){ return solventEpsStatic; }
    void setSolventEpsStatic( double epsstatic ){ solventEpsStatic = epsstatic; }
    double getSolventEpsOptical(){ return solventEpsOptical; }
    void setSolventEpsOptical( double epsoptical ){ solventEpsOptical = epsoptical; }
    double getSolventRadius(){ return solventRadius; }
    void setSolventRadius( double radius ){ solventRadius = radius; }
 private:
    string solventName;
    double solventEpsStatic;
    double solventEpsOptical;
    double solventRadius;
};

#endif
