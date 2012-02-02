/* \file GreensFunction.cpp ABC for green´s function: implemetation

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
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "MetalSphere.h"
#include "GreensFunctionSum.h"

template<class T> class Vacuum;

typedef taylor <double, 1, 1> T_DER;
typedef taylor <double, 3, 1> T_GRA;
typedef taylor <double, 3, 2> T_HES;
typedef GreensFunction<double>    G_DBL;
typedef GreensFunction<T_DER>     G_DER; 
typedef GreensFunction<T_GRA>     G_GRA;
typedef GreensFunction<T_HES>     G_HES;
typedef Vacuum<double>            V_DBL;
typedef Vacuum<T_DER>             V_DER; 
typedef Vacuum<T_GRA>             V_GRA;
typedef Vacuum<T_HES>             V_HES;
typedef UniformDielectric<double> U_DBL;
typedef UniformDielectric<T_DER>  U_DER; 
typedef UniformDielectric<T_GRA>  U_GRA;
typedef UniformDielectric<T_HES>  U_HES;
typedef GreensFunctionSum<double> S_DBL;
typedef GreensFunctionSum<T_DER>  S_DER; 
typedef GreensFunctionSum<T_GRA>  S_GRA;
typedef GreensFunctionSum<T_HES>  S_HES;

template <class T>
void GreensFunction<T>::setDelta(double value) {
	if (value <= 1.0e-10) {
		std::cout << "Delta value must be larger than 1.0e-10 " << std::endl;
		exit(-1);
	}
	delta = value;
}

template <class T>
double GreensFunction<T>::evalf(Vector3d & source, Vector3d & probe) {
	T sp[3], pp[3], res;
	sp[0] = source(0);
	sp[1] = source(1);
	sp[2] = source(2);
	pp[0] = probe(0);
	pp[1] = probe(1);
	pp[2] = probe(2);
	res = evalGreensFunction(sp, pp);
	return res[0];
}

template <>
double GreensFunction<double>::evalf(Vector3d & source, Vector3d & probe) {
	double sp[3], pp[3], res;
	sp[0] = source(0);
	sp[1] = source(1);
	sp[2] = source(2);
	pp[0] = probe(0);
	pp[1] = probe(1);
	pp[2] = probe(2);
	res = evalGreensFunction(sp, pp);
	return res;
}

template <class T>
double GreensFunction<T>::derivativeSource(Vector3d &direction, Vector3d &p1, Vector3d &p2) {
	T t1[3], t2[3], derivative;
	//	direction.normalize();
	t1[0] = p1(0); t1[0][1] = direction(0);
	t1[1] = p1(1); t1[1][1] = direction(1);
	t1[2] = p1(2); t1[2][1] = direction(2);
	t2[0] = p2(0);
	t2[1] = p2(1);
	t2[2] = p2(2);
	derivative = evalGreensFunction(t1, t2);
	return derivative[1];
}

template <>
double GreensFunction<double>::derivativeSource(Vector3d &direction, Vector3d &p1, Vector3d &p2) {
    Vector3d deltaPlus  = p1 + direction * delta / direction.norm();
    Vector3d deltaMinus = p1 - direction * delta / direction.norm();
    double funcPlus  = evalf(deltaPlus,  p2);
    double funcMinus = evalf(deltaMinus, p2);
    return (funcPlus - funcMinus)/(2.0*delta);
}

template <class T>
double GreensFunction<T>::derivativeProbe(Vector3d &direction, Vector3d &p1, Vector3d &p2) {
	T t1[3], t2[3], derivative;
	//	direction.normalize();
	t1[0] = p1(0);
	t1[1] = p1(1);
	t1[2] = p1(2);
	t2[0] = p2(0); t2[0][1] = direction(0);
	t2[1] = p2(1); t2[1][1] = direction(1);
	t2[2] = p2(2); t2[2][1] = direction(2);
	derivative = evalGreensFunction(t1, t2);
	return derivative[1];
}

template <>
double GreensFunction<double>::derivativeProbe(Vector3d &direction, Vector3d &p1, Vector3d &p2) {
    Vector3d deltaPlus  = p2 + direction * delta / direction.norm();
    Vector3d deltaMinus = p2 - direction * delta / direction.norm();
    double funcPlus  = evalf(p1, deltaPlus);
    double funcMinus = evalf(p1, deltaMinus);
    return (funcPlus - funcMinus)/(2.0*delta);
}

template <class T>
Vector3d GreensFunction<T>::gradientSource(Vector3d &p1, Vector3d &p2) {
	Vector3d g;
	gradientSource(g, p1, p2);
    return g;
}

template <class T>
void GreensFunction<T>::gradientSource(Vector3d &g, Vector3d &p1, Vector3d &p2) {
	T t1[3], t2[3], grad;
	t1[0] = p1(0); t1[0][1] = 1;
	t1[1] = p1(1); t1[1][2] = 1;
	t1[2] = p1(2); t1[2][3] = 1;
	t2[0] = p2(0);
	t2[1] = p2(1);
	t2[2] = p2(2);
	grad = evalGreensFunction(t1, t2);
    g << grad[1], grad[2], grad[3];
    return;
}

template <>
void GreensFunction<double>::gradientSource(Vector3d &g, Vector3d &p1, Vector3d &p2) {
	Vector3d direction;
	direction << 1.0, 0.0, 0.0;
	g(0) = derivativeSource(direction, p1, p2);
	direction << 0.0, 1.0, 0.0;
	g(1) = derivativeSource(direction, p1, p2);
	direction << 0.0, 0.0, 1.0;
	g(2) = derivativeSource(direction, p1, p2);
}

template <>
void GreensFunction< taylor<double, 1, 1> >::gradientSource(Vector3d &g, Vector3d &p1, Vector3d &p2) {
	Vector3d direction;
	direction << 1.0, 0.0, 0.0;
	g(0) = derivativeSource(direction, p1, p2);
	direction << 0.0, 1.0, 0.0;
	g(1) = derivativeSource(direction, p1, p2);
	direction << 0.0, 0.0, 1.0;
	g(2) = derivativeSource(direction, p1, p2);
}

template <class T>
Vector3d GreensFunction<T>::gradientProbe(Vector3d &p1, Vector3d &p2) {
	Vector3d g;
	gradientProbe(g, p1, p2);
    return g;
}

template <class T>
void GreensFunction<T>::gradientProbe(Vector3d &g, Vector3d &p1, Vector3d &p2) {
	T t1[3], t2[3], grad;
	t1[0] = p1(0);
	t1[1] = p1(1);
	t1[2] = p1(2);
	t2[0] = p2(0); t2[0][1] = 1;
	t2[1] = p2(1); t2[1][2] = 1;
	t2[2] = p2(2); t2[2][3] = 1;
	grad = evalGreensFunction(t1, t2);
    g << grad[1], grad[2], grad[3];
    return;
}

template <>
void GreensFunction<double>::gradientProbe(Vector3d &g, Vector3d &p1, Vector3d &p2) {
	Vector3d direction;
	direction << 1.0, 0.0, 0.0;
	g(0) = derivativeProbe(direction, p1, p2);
	direction << 0.0, 1.0, 0.0;
	g(1) = derivativeProbe(direction, p1, p2);
	direction << 0.0, 0.0, 1.0;
	g(2) = derivativeProbe(direction, p1, p2);
}

template <>
void GreensFunction<taylor<double, 1, 1> >::gradientProbe(Vector3d &g, Vector3d &p1, Vector3d &p2) {
	Vector3d direction;
	direction << 1.0, 0.0, 0.0;
	g(0) = derivativeProbe(direction, p1, p2);
	direction << 0.0, 1.0, 0.0;
	g(1) = derivativeProbe(direction, p1, p2);
	direction << 0.0, 0.0, 1.0;
	g(2) = derivativeProbe(direction, p1, p2);
}

template <class T>
GreensFunction<T>* GreensFunction<T>::allocateGreensFunction(const Section &green) {
	GreensFunction<T> *gf = 0;
	const string greenType = green.getStr("Type");
	if (greenType == "Vacuum") {
		gf = new Vacuum<T>();
	} else if (greenType == "UniformDielectric") {
		gf = new UniformDielectric<T>(green);
	} else if (greenType == "GreensFunctionSum") {
		gf = new GreensFunctionSum<T>(green);
	} else {
		cout << "Unknown Greens function" << endl;
		exit(1);
	}
	return gf;
}

template <>
GreensFunction<double>* GreensFunction<double>::allocateGreensFunction(const Section &green) {
	GreensFunction<double> *gf;
	const string greenType = green.getStr("Type");
	if (greenType == "Vacuum") {
		gf = new Vacuum<double>();
	} else if (greenType == "UniformDielectric") {
		gf = new UniformDielectric<double>(green);
	} else if (greenType == "MetalSphere") {
		gf = new MetalSphere(green);
	} else if (greenType == "GreensFunctionSum") {
		gf = new GreensFunctionSum<double>(green);
	} else {
		cout << "Unknown Greens function" << endl;
		exit(1);
	}
	return gf;
}

template <class T>
GreensFunction<T>* GreensFunction<T>::allocateGreensFunction(double dielConst) {
	GreensFunction<T> *gf;
	gf = new UniformDielectric<T>(dielConst);
	return gf;
}

template <class T>
GreensFunction<T>* GreensFunction<T>::allocateGreensFunction() {
	GreensFunction<T> *gf;
	gf = new Vacuum<T>();
	return gf;
}

template <class T>
std::ostream & GreensFunction<T>::printObject(std::ostream &os) {
	os << "Green's Function" << std::endl;
	os << "Delta = " << delta << std::endl;
	os << "Uniform = " << uniformFlag; 
	return os;
}

template class GreensFunction<double>;
template class GreensFunction<taylor <double, 1, 1> >;
template class GreensFunction<taylor <double, 3, 1> >;
template class GreensFunction<taylor <double, 3, 2> >;
template class Vacuum<double>;
template class Vacuum<taylor <double, 1, 1> >;
template class Vacuum<taylor <double, 3, 1> >;
template class Vacuum<taylor <double, 3, 2> >;
