#include <iostream>
#include <ostream>
#include <cmath>
#include <cstdlib>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "MetalSphere.h"
#include "GreensFunctionSum.h"

GreensFunctionSum::GreensFunctionSum(GreensFunction &first, 
									 GreensFunction &second){
	greenFirst  = &first;
	greenSecond = &second;
	uniformFlag = greenFirst->isUniform() && greenSecond->isUniform();
};

GreensFunctionSum::GreensFunctionSum(Section green){
	greenFirst  = allocateGreensFunction(green.getSect("Green<one>"));
	greenSecond = allocateGreensFunction(green.getSect("Green<two>"));
	uniformFlag = greenFirst->isUniform() && greenSecond->isUniform();
};

double GreensFunctionSum::evalf(Vector3d &p1, Vector3d &p2) {
    double valFirst = greenFirst->evalf(p1, p2);
    double valSecond = greenSecond->evalf(p1, p2);
    return valFirst + valSecond;
}

double GreensFunctionSum::derivative(Vector3d &direction, Vector3d &p1, Vector3d &p2) {
    double derFirst = greenFirst->derivative(direction, p1, p2);
    double derSecond = greenSecond->derivative(direction, p1, p2);
    return derFirst + derSecond;
}

double GreensFunctionSum::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2) {
    double derFirst = greenFirst->evald(direction, p1, p2);
    double derSecond = greenSecond->evald(direction, p1, p2);
    return derFirst + derSecond;
}

void GreensFunctionSum::gradient(Vector3d &gradient, Vector3d &p1, Vector3d &p2) {
    Vector3d gradFirst, gradSecond;
    greenFirst->gradient(gradFirst, p1, p2);
    greenSecond->gradient(gradSecond, p1, p2);
    gradient = gradFirst + gradSecond;
    return;
}

