/* \file GreensFunction.cpp ABC for green´s function: implemetation

 */
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

#include "GreensFunction.h"

/** Computes numerically the directional derivative of the Green's function

    @param p1 source point
    @param p2 potential point
    @param delta derivative step

 */
double GreensFunction::derivative(Vector3d &direction, Vector3d &p1, Vector3d &p2, double delta) {
    Vector3d deltaPlus, deltaMinus;
    double norm = sqrt(direction.dot(direction));
    deltaPlus  = p1 + direction * delta / norm;
    deltaMinus = p1 - direction * delta / norm;
    double funcPlus  = evalf(deltaPlus,  p2);
    double funcMinus = evalf(deltaMinus, p2);
    double numDer = (funcPlus - funcMinus)/(2.0*delta);
    return numDer;
}

void GreensFunction::gradient(Vector3d &gradient, Vector3d &p1, Vector3d &p2, double delta) {
    Vector3d xdir(1.0, 0.0, 0.0);
    Vector3d ydir(0.0, 1.0, 0.0);
    Vector3d zdir(0.0, 0.0, 1.0);
    gradient(0) = derivative(xdir, p1, p2, delta);
    gradient(1) = derivative(ydir, p1, p2, delta);
    gradient(2) = derivative(zdir, p1, p2, delta);
    return;
}

