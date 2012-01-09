#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Solvent.h"

/*

  Methods for Solvent class
  written by Roberto Di Remigio, 2011

*/

Solvent::Solvent( const string & name, double epsstatic, 
				  double epsoptical, double radius ) {
	solventName = name;
	solventEpsStatic = epsstatic;
	solventEpsOptical = epsoptical;
	solventRadius = radius;
}
