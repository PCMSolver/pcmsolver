/* \file GreensFunction.cpp ABC for green´s function: implemetation

 */
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

#include "Getkw.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "MetalSphere.h"
#include "GreensFunctionSum.h"


/** Computes numerically the directional derivative of the Green's function

    @param p1 source point
    @param p2 potential point
    @param delta derivative step

 */

void GreensFunction::setDelta(double value) {
	if (value <= 1.0e-10) {
		cout << "Delta value must be larger than 1.0e-10 " << endl;
		exit(-1);
	}
	delta = value;
}

double GreensFunction::derivative(Vector3d &direction, Vector3d &p1, Vector3d &p2) {
    Vector3d deltaPlus, deltaMinus;
    double norm = sqrt(direction.dot(direction));
    deltaPlus  = p1 + direction * delta / norm;
    deltaMinus = p1 - direction * delta / norm;
    double funcPlus  = evalf(deltaPlus,  p2);
    double funcMinus = evalf(deltaMinus, p2);
    double numDer = (funcPlus - funcMinus)/(2.0*delta);
    return numDer;
}

void GreensFunction::gradient(Vector3d &gradient, Vector3d &p1, Vector3d &p2) {
    Vector3d xdir(1.0, 0.0, 0.0);
    Vector3d ydir(0.0, 1.0, 0.0);
    Vector3d zdir(0.0, 0.0, 1.0);
    gradient(0) = derivative(xdir, p1, p2);
    gradient(1) = derivative(ydir, p1, p2);
    gradient(2) = derivative(zdir, p1, p2);
    return;
}

GreensFunction* GreensFunction::allocateGreensFunction(const Section &green) {
	GreensFunction *gf;
	const string greenType = green.getStr("Type");
	if (greenType == "Vacuum") {
		gf = new Vacuum();
	} else if (greenType == "UniformDielectric") {
		gf = new UniformDielectric(green);
	} else if (greenType == "MetalSphere") {
		gf = new MetalSphere(green);
	} else if (greenType == "GreensFunctionSum") {
		gf = new GreensFunctionSum(green);
	} else {
		cout << "Unknown Greens function" << endl;
		exit(1);
	}
	return gf;
}

