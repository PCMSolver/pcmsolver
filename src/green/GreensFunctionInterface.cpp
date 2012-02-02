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
#include "GreensFunction.h"


typedef taylor <double, 1, 1> T_DERIVATIVE;
typedef taylor <double, 3, 1> T_GRADIENT;
typedef taylor <double, 3, 2> T_HESSIAN;
typedef GreensFunction<double> G_DOUBLE;
typedef GreensFunction<T_DERIVATIVE>  G_DERIVATIVE; 
typedef GreensFunction<T_GRADIENT>  G_GRADIENT;
typedef GreensFunction<T_HESSIAN>  G_HESSIAN;

G_DOUBLE     * g0;
G_DERIVATIVE * g1;
G_GRADIENT   * g2;
G_HESSIAN    * g3;

GreensFunctionInterface* 
GreensFunctionInterface::allocateGreensFunctionInterface(const Section &green) {
	GreensFunctionInterface *gf;
	const string greenDer = green.getStr("Der");
	if (greenDer == "Numerical") {
		gf = g0->allocateGreensFunction(green);
	} else if (greenDer == "Derivative") {
		gf = g1->allocateGreensFunction(green);
	} else if (greenDer == "Gradient") {
		gf = g2->allocateGreensFunction(green);
	} else if (greenDer == "Hessian") {
		gf = g3->allocateGreensFunction(green);
	}
	return gf;
}

GreensFunctionInterface* 
GreensFunctionInterface::allocateGreensFunctionInterface(double epsilon, const string greenDer) {
 	GreensFunctionInterface *gf;
	if (greenDer == "Numerical") {
		gf = g0->allocateGreensFunction(epsilon);
	} else if (greenDer == "Derivative") {
		gf = g1->allocateGreensFunction(epsilon);
	} else if (greenDer == "Gradient") {
		gf = g2->allocateGreensFunction(epsilon);
	} else if (greenDer == "Hessian") {
		gf = g3->allocateGreensFunction(epsilon);
	}
	return gf;
}

GreensFunctionInterface* 
GreensFunctionInterface::allocateGreensFunctionInterface(const string greenDer) {
	GreensFunctionInterface *gf;
	if (greenDer == "Numerical") {
		gf = g0->allocateGreensFunction();
	} else if (greenDer == "Derivative") {
		gf = g1->allocateGreensFunction();
	} else if (greenDer == "Gradient") {
		gf = g2->allocateGreensFunction();
	} else if (greenDer == "Hessian") {
		gf = g3->allocateGreensFunction();
	}
	return gf;
}

/*
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
*/
