#include <iostream>
#include <ostream>
#include <cmath>
#include <cstdlib>
#include <Eigen/Dense>

using namespace Eigen;

using namespace std;

#include "GreensFunction.h"
#include "MetalSphere.h"

extern"C" {
    void gsfera_cpp_(double* epssol, double* epsre, double* epsim, 
		     double* sphRadius, double* ps, double* p1, double* p2,
		     double* greenre, double* greenim);
}

MetalSphere::MetalSphere(double eps, double epsRe, double epsIm, Vector3d &pos,
			 double radius){
        epsSolvent = eps;
        sphRadius = radius;
        sphPosition = pos;
        epsMetal = dcomplex(epsRe,epsIm);
        uniformFlag = false;
    };

double MetalSphere::evalf(Vector3d &p1, Vector3d &p2) {
    double epsre, epsim, greenre, greenim;
    double point1[3], point2[3], sphere[3];
    point1[0] = p1(0);
    point1[1] = p1(1);
    point1[2] = p1(2);
    point2[0] = p2(0);
    point2[1] = p2(1);
    point2[2] = p2(2);
    sphere[0] = sphPosition(0);
    sphere[1] = sphPosition(1);
    sphere[2] = sphPosition(2);
    epsre =  epsMetal.real();
    epsim =  epsMetal.imag();
    gsfera_cpp_(&epsSolvent, &epsre, &epsim, &sphRadius,
	       sphere, point1, point2, &greenre, &greenim);
    return greenre;
}

double MetalSphere::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2, double delta) {
    return epsSolvent * derivative(direction, p1, p2, delta);
}

