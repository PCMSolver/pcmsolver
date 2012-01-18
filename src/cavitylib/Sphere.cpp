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

Sphere::Sphere( double radius, Vector3d & center) {
  sphereRadius = radius;
  sphereCenter = center;
}

Sphere::Sphere( Atom & atom ) {
  sphereRadius = atom.getAtomRadius();
  sphereCenter = atom.getAtomCoord();
}
