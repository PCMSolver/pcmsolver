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

G_DOUBLE     * g0 = 0;
G_DERIVATIVE * g1 = 0;
G_GRADIENT   * g2 = 0;
G_HESSIAN    * g3 = 0;

GreensFunctionInterface* 
GreensFunctionInterface::allocateGreensFunctionInterface(const Section &green) 
{
	std::cout << "Allocator from input section" << std::endl;
	GreensFunctionInterface *gf = 0;
	const string greenDer = green.getStr("Der");
	if (greenDer == "Numerical") {
		gf = g0->allocateGreensFunction(green);
		std::cout << "Numerical der" << std::endl;
	} else if (greenDer == "Derivative") {
		gf = g1->allocateGreensFunction(green);
		std::cout << "Analytical der" << std::endl;
	} else if (greenDer == "Gradient") {
		gf = g2->allocateGreensFunction(green);
		std::cout << "Full gradient der" << std::endl;
	} else if (greenDer == "Hessian") {
		gf = g3->allocateGreensFunction(green);
		std::cout << "Full hessian der" << std::endl;
	}
	return gf;
}

GreensFunctionInterface* 
GreensFunctionInterface::allocateGreensFunctionInterface(double epsilon, const std::string greenDer) 
{
	std::cout << "Allocator from epsilon and der type." << std::endl;
 	GreensFunctionInterface *gf = 0;
	if (greenDer == "Numerical") {
		gf = g0->allocateGreensFunction(epsilon);
		std::cout << "Numerical der" << std::endl;
	} else if (greenDer == "Derivative") {
		gf = g1->allocateGreensFunction(epsilon);
		std::cout << "Analytical der" << std::endl;
	} else if (greenDer == "Gradient") {
		gf = g2->allocateGreensFunction(epsilon);
		std::cout << "Full gradient der" << std::endl;
	} else if (greenDer == "Hessian") {
		gf = g3->allocateGreensFunction(epsilon);
		std::cout << "Full hessian der" << std::endl;
	}
	return gf;
}

GreensFunctionInterface* 
GreensFunctionInterface::allocateGreensFunctionInterface(const std::string greenDer) 
{
	std::cout << "Allocator from der type." << std::endl;
	GreensFunctionInterface *gf = 0;
	if (greenDer == "Numerical") {
		gf = g0->allocateGreensFunction();
		std::cout << "Numerical der" << std::endl;
	} else if (greenDer == "Derivative") {
		gf = g1->allocateGreensFunction();
		std::cout << "Analytical der" << std::endl;
	} else if (greenDer == "Gradient") {
		gf = g2->allocateGreensFunction();
		std::cout << "Full gradient der" << std::endl;
	} else if (greenDer == "Hessian") {
		gf = g3->allocateGreensFunction();
		std::cout << "Full hessian der" << std::endl;
	}
	return gf;
}
