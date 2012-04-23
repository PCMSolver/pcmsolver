#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Atom.h"

/*

  Methods for Atom class
  written by Roberto Di Remigio, 2011

*/

Atom::Atom( const string & element, const string & symbol, double charge, 
			double radius, Vector3d & coord, double scaling, const string & colour ) {
  atomElement = element;
  atomSymbol = symbol;
  atomCharge = charge;
  atomRadius = radius;
  atomCoord = coord;
  atomColour = colour;
  atomRadiusScaling = scaling;
}

Atom::Atom( const string & element, const string & symbol, double charge, 
			double radius ) {
  Vector3d Origin(0.0, 0.0, 0.0);
  string colour = "Violet";
  atomElement = element;
  atomSymbol = symbol;
  atomCharge = charge;
  atomRadius = radius;
  atomCoord = Origin;
  atomColour = colour;
  atomRadiusScaling = 1.0;
}
