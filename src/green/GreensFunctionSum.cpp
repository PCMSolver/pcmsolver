#include <iostream>
#include <ostream>
#include <cmath>
#include <cstdlib>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "GreensFunction.h"
#include "GreensFunctionSum.h"

GreensFunctionSum::GreensFunctionSum(GreensFunction &first, GreensFunction &second){
        greenFirst  = &first;
        greenSecond = &second;
	uniformFlag = greenFirst->isUniform() && greenSecond->isUniform();
    };

double GreensFunctionSum::evalf(Vector3d &p1, Vector3d &p2) {
    double valFirst = greenFirst->evalf(p1, p2);
    double valSecond = greenSecond->evalf(p1, p2);
    return valFirst + valSecond;
}

double GreensFunctionSum::derivative(Vector3d &direction, Vector3d &p1, Vector3d &p2, double delta) {
    double derFirst = greenFirst->derivative(direction, p1, p2, delta);
    double derSecond = greenSecond->derivative(direction, p1, p2, delta);
    return derFirst + derSecond;
}

double GreensFunctionSum::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2, double delta) {
    double derFirst = greenFirst->evald(direction, p1, p2, delta);
    double derSecond = greenSecond->evald(direction, p1, p2, delta);
    return derFirst + derSecond;
}

void GreensFunctionSum::gradient(Vector3d &gradient, Vector3d &p1, Vector3d &p2, double delta) {
    Vector3d gradFirst, gradSecond;
    greenFirst->gradient(gradFirst, p1, p2, delta);
    greenSecond->gradient(gradSecond, p1, p2, delta);
    gradient = gradFirst + gradSecond;
    return;
}
