#ifndef SPHERE_H
#define SPHERE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Config.h"

#include <Eigen/Dense>

#include "Atom.h"

using namespace std;
using namespace Eigen;

/*

  A Sphere class
  written by Roberto Di Remigio, 2011

*/

class Sphere
{
 private:
  double sphereRadius;
  Vector3d sphereCenter;
  string sphereColour;
       
 protected:
  virtual ostream & printObject(ostream & os);

 public:
  Sphere(){}
  Sphere( Vector3d & sphereCenter, double sphereRadius, const string & sphereColour = "Violet" );
  Sphere( Atom & atom );
  //  Sphere( Atom & atom, double charge );
  ~Sphere(){}
  double getSphereRadius(){ return sphereRadius; }
  void setSphereRadius( double radius ){ sphereRadius = radius; }
  Vector3d & getSphereCenter(){ return sphereCenter; }
  double getSphereCenter(int i){ return sphereCenter(i); }
  void setSphereCenter( Vector3d & coord ){ sphereCenter = coord; }
  string & getSphereColour(){ return sphereColour; }
  void setSphereColour( string & colour ){ sphereColour = colour; }
  
  friend std::ostream& operator<<(std::ostream & o, Sphere & s);
};

#endif
