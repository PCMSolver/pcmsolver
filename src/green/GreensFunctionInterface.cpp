/* \file GreensFunctionInterface.cpp ABC for green´s function: implemetation

 */
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

#include "Getkw.h"
#include "taylor.hpp"
//#include "TaylorSupport.h"
#include "GreensFunctionInterface.h"


template class taylor <double, 1, 1>;
template class taylor <double, 3, 1>;
template class taylor <double, 3, 2>;
template class GreensFunction<double>;
template class GreensFunction<taylor <double, 1, 1> >;
template class GreensFunction<taylor <double, 3, 1> >;
template class GreensFunction<taylor <double, 3, 2> >;

GreensFunctionInterface* 
GreensFunctionInterface::allocateGreensFunctionInterface(const Section &green) {
	GreensFunctionInterface *gfi;
	const string greenDer = green.getStr("Der");
	if (greenDer == "Numerical") {
		gf = allocateGreensFunction<double>(const Section &green);
	} else if (greenDer == "Derivative") {
		gf = allocateGreensFunction< taylor <double, 1, 1> >(const Section &green);
	} else if (greenDer == "Gradient") {
		gf = allocateGreensFunction< taylor <double, 3, 1> >(const Section &green);
	} else if (greenDer == "Hessian") {
		gf = allocateGreensFunction< taylor <double, 3, 2> >(const Section &green);
	}
	return gf;
}

GreensFunctionInterface* 
GreensFunctionInterface::allocateGreensFunctionInterface(const string greenDer="Derivative", double epsilon) { 
	GreensFunctionInterface *gfi;
	const string greenDer = green.getStr("Derivative");
	if (greenDer == "Numerical") {
		gf = allocateGreensFunction<double>(epsilon);
	} else if (greenDer == "Derivative") {
		gf = allocateGreensFunction< taylor <double, 1, 1> >(epsilon);
	} else if (greenDer == "Gradient") {
		gf = allocateGreensFunction< taylor <double, 3, 1> >(epsilon);
	} else if (greenDer == "Hessian") {
		gf = allocateGreensFunction< taylor <double, 3, 2> >(epsilon);
	}
	return gf;
}

GreensFunctionInterface* 
GreensFunctionInterface::allocateGreensFunctionInterface(const string greenDer="Derivative") {
	GreensFunctionInterface *gfi;
	if (greenDer == "Numerical") {
		gf = allocateGreensFunction<double>();
	} else if (greenDer == "Derivative") {
		gf = allocateGreensFunction< taylor <double, 1, 1> >();
	} else if (greenDer == "Gradient") {
		gf = allocateGreensFunction< taylor <double, 3, 1> >();
	} else if (greenDer == "Hessian") {
		gf = allocateGreensFunction< taylor <double, 3, 2> >();
	}
	return gf;
}

