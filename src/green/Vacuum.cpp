#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

#include "taylor.hpp"
//#include "TaylorSupport.h"
#include "GreensFunction.h"
#include "Vacuum.h"

class Section;

template<class T>
T Vacuum<T>::evalGreensFunction(T * sp, T * pp) {
	T res;
	res = 1.0/sqrt((sp[0]-pp[0])*(sp[0]-pp[0])+
				   (sp[1]-pp[1])*(sp[1]-pp[1])+
				   (sp[2]-pp[2])*(sp[2]-pp[2]));
	return res;
}

template<class T>
double Vacuum<T>::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2){
    return direction.dot(this->gradientProbe(p1, p2))/direction.norm();
}

template class Vacuum<double>;
template class Vacuum< taylor<double, 3, 1> >;
template class Vacuum< taylor<double, 3, 2> >;
