#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

#include "GreensFunction.h"
#include "Vacuum.h"

double Vacuum::evalf(Vector3d &p1, Vector3d &p2) {
    double dist = sqrt((p1 - p2).dot(p1 - p2));
    return 1.0 / dist;
}

void Vacuum::gradient(Vector3d &gradient, Vector3d &p1, Vector3d &p2, double delta){
    double dist = sqrt((p1 - p2).dot(p1 - p2));
    double distcube = dist * dist * dist;
    gradient = (p2 - p1) / distcube;
    std::cout << "analytical gradient" << std::endl;
    return;
}

double Vacuum::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2, double delta){
    double dist = sqrt((p1 - p2).dot(p1 - p2));
    double distcube = dist * dist * dist;
    double der = direction.dot(p2 - p1) / distcube;
    return der;
}
double Vacuum::derivative(Vector3d &direction, Vector3d &p1, Vector3d &p2, double delta){
    double dist = sqrt((p1 - p2).dot(p1 - p2));
    double distcube = dist * dist * dist;
    double der = direction.dot(p2 - p1) / distcube;
    std::cout << "analytical der" << std::endl;
    return der;
}

