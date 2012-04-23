#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

#include "SurfaceFunction.h"

using namespace std;
using namespace Eigen;


SurfaceFunction::SurfaceFunction(const std::string & name) {
	this->name = name;
	allocated = false;
}

SurfaceFunction::SurfaceFunction(const std::string & name, int nPoints) {
	this->name = name;
	this->values.resize(nPoints);
	allocated = true;
}

SurfaceFunction::SurfaceFunction(const std::string & name, int nPoints, double * values) {
	this->name = name;
	this->values.resize(nPoints);
	allocated = true;
	for (int i = 0; i < nPoints; i++) {
		this->values(i) = values[i];
	}
}

void SurfaceFunction::setValues(double * values) {
	for (int i = 0; i < this->values.size(); i++) {
		this->values(i) = values[i];
	}
}

void SurfaceFunction::getValues(double * values) {
	for (int i = 0; i < this->values.size(); i++) {
		values[i] = this->values(i);
	}
}

void SurfaceFunction::clear() {
	this->values.setZero();
}

ostream & operator<<(ostream & os, SurfaceFunction & sf) {
	return sf.printObject(os);
}

ostream & SurfaceFunction::printObject(ostream & os) {
	os << "Surface Function " << this->name << endl;
	if (this->isAllocated()) {
		os << values.transpose() << endl;
	}
	return os;
}

