#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

#include "Getkw.h"
#include "GreensFunction.h"
#include "UniformDielectric.h"

UniformDielectric::UniformDielectric(double dielConst) 
{
    epsilon = dielConst;
    uniformFlag = true;
}

UniformDielectric::UniformDielectric(Section green) 
{
    epsilon = green.getDbl("Eps");
    uniformFlag = true;
}

double UniformDielectric::evalf(Vector3d &p1, Vector3d &p2) {
    double dist = sqrt((p1 - p2).dot(p1 - p2));
    return 1.0/(epsilon * dist);
}

void UniformDielectric::gradient(Vector3d &gradient, Vector3d &p1, Vector3d &p2){
    double dist = sqrt((p1 - p2).dot(p1 - p2));
    double distcube = dist * dist * dist;
    gradient = (p2 - p1) / (epsilon * distcube);
    std::cout << "analytical gradient" << std::endl;
    return;
}

double UniformDielectric::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2){
    double dist = sqrt((p1 - p2).dot(p1 - p2));
    double distcube = dist * dist * dist;
    double der = direction.dot(p2 - p1) / distcube;
    return der;
}
double UniformDielectric::derivative(Vector3d &direction, Vector3d &p1, Vector3d &p2){
    double dist = sqrt((p1 - p2).dot(p1 - p2));
    double distcube = dist * dist * dist;
    double der = direction.dot(p2 - p1) / (epsilon * distcube);
    std::cout << "analytical der" << std::endl;
    return der;
}

