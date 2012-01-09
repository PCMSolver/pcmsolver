#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Atom.h"
#include "Sphere.h"

/*

  Methods for Sphere class
  written by Roberto Di Remigio, 2011

*/

Sphere::Sphere( Vector3d & center, double radius ) {
	sphereCenter = center;
	sphereRadius = radius;
}

Sphere::Sphere( Atom & atom ) {
	sphereCenter = atom.getAtomCoord();
	sphereRadius = atom.getAtomRadius();
}

ostream & operator<<(ostream & os, Sphere & sphere) {
	return sphere.printObject(os);
} 

ostream & Sphere::printObject(ostream & os) {
	os << "Sphere radius " << sphereRadius << endl;
	os << "Sphere center\n" << sphereCenter;
	return os;
}
