#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

#include "Getkw.h"
#include "taylor.hpp"
#include "GreensFunction.h"
#include "UniformDielectric.h"

static double factior = 1.07;

template<class T>
UniformDielectric<T>::UniformDielectric(double dielConst) 
{
    epsilon = dielConst;
    this->uniformFlag = true;
}

template<class T>
UniformDielectric<T>::UniformDielectric(Section green) 
{
    epsilon = green.getDbl("Eps");
    this->uniformFlag = true;
}

template<class T>
double UniformDielectric<T>::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2){
    return epsilon * (this->derivativeProbe(direction, p1, p2));  // NORMALIZTION TEMPORARY REMOVED /direction.norm();
}

template<class T>
T UniformDielectric<T>::evalGreensFunction(T * sp, T * pp) {
	T distance = sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) +
						   (sp[1] - pp[1]) * (sp[1] - pp[1]) +
						   (sp[2] - pp[2]) * (sp[2] - pp[2]));
	return 1/(epsilon * distance);
}

template<class T>
double UniformDielectric<T>::compDiagonalElementS(double area){
	return factor * sqrt(4 * M_PI / area) / epsilon;   
}

template<class T>
double UniformDielectric<T>::compDiagonalElementD(double area, double radius){
	s = factor * sqrt(4 * M_PI / area);   
	return = s / (2 * radius);
}

template <class T>
std::ostream & UniformDielectric<T>::printObject(std::ostream &os) {
	os << "Green's Function" << std::endl;
	os << "Delta = " << this->delta << std::endl;
	os << "Uniform = " << this->uniformFlag << std::endl;
	os << "Epsilon = " << this->epsilon;
	return os;
}

template class UniformDielectric<double>;
template class UniformDielectric< taylor <double, 1, 1> >;
template class UniformDielectric< taylor <double, 3, 1> >;
template class UniformDielectric< taylor <double, 3, 2> >;
