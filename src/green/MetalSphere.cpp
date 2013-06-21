#include <iostream>
#include <ostream>
#include <cmath>
#include <cstdlib>
#include <Eigen/Dense>

using namespace Eigen;

using namespace std;

#include "Getkw.h"
#include "GreensFunctionInterface.h"
#include "GreensFunction.h"
#include "MetalSphere.h"

extern"C" {
    void gsfera_cpp_(double* epssol, double* epsre, double* epsim, 
		     double* sphRadius, double* ps, double* p1, double* p2,
		     double* greenre, double* greenim);
}

/*
MetalSphere::MetalSphere(Section green){
	epsSolvent = green.getDbl("Eps");
	epsMetal = dcomplex(green.getDbl("EpsRe"), green.getDbl("EpsImg"));
	sphRadius = green.getDbl("SphereRadius");
	const vector<double> &pos_ = green.getDblVec("SpherePosition");
	Vector3d pos(pos_[0], pos_[1], pos_[2]);
	uniformFlag = false;
};*/

double MetalSphere::evalGreensFunction(double * source, double * probe) {
    double epsre, epsim;
	double greenre, greenim;
    double point1[3], point2[3], sphere[3];
    point1[0] = source[0];
    point1[1] = source[1];
    point1[2] = source[2];
    point2[0] =  probe[0];
    point2[1] =  probe[1];
    point2[2] =  probe[2];
    sphere[0] = sphPosition(0);
    sphere[1] = sphPosition(1);
    sphere[2] = sphPosition(2);
    epsre = epsMetal.real();
    epsim = epsMetal.imag();
    gsfera_cpp_(&epsSolvent, &epsre, &epsim, &sphRadius,
	       sphere, point1, point2, &greenre, &greenim);
    return greenre;
}

double MetalSphere::compDiagonalElementS(double area){
	std::cout << "Not Yet Implemented" << std::endl;
	exit(-1);
	return 0;
}

double MetalSphere::compDiagonalElementD(double area, double radius){
	std::cout << "Not Yet Implemented" << std::endl;
	exit(-1);
	return 0;
}

double MetalSphere::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2) {
    return epsSolvent * (this->derivativeProbe(direction, p1, p2));  // NORMALIZTION TEMPORARY REMOVED /direction.norm();
}

