#ifndef ATOM_H
#define ATOM_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "Constants.h"

using namespace std;
using namespace Eigen;

/*

  An Atom class
  written by Roberto Di Remigio, 2011

  "Why use a structure when you can define a class?"

  vector<Atom> Bondi[] should contain the van der Waals radii taken from 
  A. Bondi, J. Phys. Chem. 68, 441-451 (1964).
  vector<Atom> Tatewaki[] should contain the effective atomic radii taken
  from H. Tatewaki et al. Bull. Chem. Soc. Jpn. 83, 1203-1210 (2010).

  They should be declared as constant vectors, with all the atomCoord 
  data members set to void.

*/

class Atom {
 private:
  string atomElement;
  string atomSymbol;
  Vector3d atomCoord;
  double atomCharge;
  double atomRadius;
  string atomColour;
  double atomRadiusScaling;
       
 public:
  Atom(){}
  Atom( const string & element, const string & symbol, double charge, 
        double radius, Vector3d & coord, double scaling = 1.0, const string & colour = "Violet" );
  Atom( const string & element, const string & symbol, double charge, 
        double radius );
  ~Atom(){}
  string getAtomElement(){ return atomElement; }
  void setAtomElement( const string & element ){ atomElement = element; }
  string getAtomSymbol(){ return atomSymbol; }
  void setAtomSymbol( const string & symbol ){ atomSymbol = symbol; }
  Vector3d getAtomCoord(){ return atomCoord; }
  void setAtomCoord( Vector3d & coord ){ atomCoord = coord; }
  double getAtomCharge(){ return atomCharge; }
  void setAtomCharge( double charge ){ atomCharge = charge; }
  double getAtomRadius(){ return (atomRadius * ToAngstrom); }
  void setAtomRadius( double radius ){ atomRadius = radius; }
  double getAtomRadiusScaling(){ return atomRadiusScaling; }
  void setAtomRadiusScaling( double scaling ){ atomRadiusScaling = scaling; }
  string getAtomColour(){ return atomColour; }
  void setAtomColour( const string & colour ){ atomColour = colour; }
  static vector<Atom> initBondi();
};

#endif
