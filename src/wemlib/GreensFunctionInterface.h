#ifndef GREENSFUNCTIONINTERFACE_H
#define GREENSFUNCTIONINTERFACE_H

// Interface for the Green's function

#include <Eigen/Core>

using namespace Eigen;

class GreensFunctionInterface{
public:

  virtual ~GreensFunctionInterface() {};

  virtual double getSingleLayer(Vector3d &x, Vector3d &y) = 0;
  virtual double getDoubleLayer(Vector3d &x, Vector3d &y, Vector3d &n_y) = 0;

};

#endif

