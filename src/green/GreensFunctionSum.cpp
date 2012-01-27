#include <iostream>
#include <ostream>
#include <cmath>
#include <cstdlib>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "taylor.hpp"
#include "GreensFunction.h"
#include "GreensFunctionSum.h"

template<class T>
GreensFunctionSum<T>::GreensFunctionSum(GreensFunction<T> &first, 
										GreensFunction<T> &second)
{
	this->greenFirst  = &first;
	this->greenSecond = &second;
	this->uniformFlag = this->greenFirst->isUniform() && 
		                this->greenSecond->isUniform();
}

template<class T>
GreensFunctionSum<T>::GreensFunctionSum(Section green){
	greenFirst  = this->allocateGreensFunction(green.getSect("Green<one>"));
	greenSecond = this->allocateGreensFunction(green.getSect("Green<two>"));
	this->uniformFlag = greenFirst->isUniform() && greenSecond->isUniform();
}

template<class T>
double GreensFunctionSum<T>::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2) {
    double derFirst = greenFirst->evald(direction, p1, p2);
    double derSecond = greenSecond->evald(direction, p1, p2);
    return derFirst + derSecond;
}

template<class T>
T GreensFunctionSum<T>::evalGreensFunction(T * sp, T * pp) {
	std::cout << "NYI" << std::endl;
	exit(-1);
	//	return greenFirst->evalGreensFunction(sp, pp)
	//		+ greenSecond->evalGreensFunction(sp, pp);
}

template class GreensFunctionSum<double>;
template class GreensFunctionSum< taylor <double, 1, 1> >;
template class GreensFunctionSum< taylor <double, 3, 1> >;
template class GreensFunctionSum< taylor <double, 3, 2> >;
