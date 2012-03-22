#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

#include "SurfaceFunction.h"

using namespace std;
using namespace Eigen;


SurfaceFunction::SurfaceFunction(std::string & name) {
	this->name = name;
	allocated = false;
}

SurfaceFunction::SurfaceFunction(std::string & name, int nPoints) {
	this->name = name;
	this->values.resize(nPoints);
	allocated = true;
}

SurfaceFunction::SurfaceFunction(std::string & name, int nPoints, double * values) {
	this->name = name;
	this->values.resize(nPoints);
	allocated = true;
	for (int i = 0; i < nPoints; i++) {
		this->values(i) = values[i];
	}
}

void SurfaceFunction::setValues(double * values) {
	for (int i = 0; i < functionValues.size(); i++) {
		this->values(i) = values[i];
	}
}

void SurfaceFunction::getValues(double * values) {
	for (int i = 0; i < functionValues.size(); i++) {
		values[i] = functionValues(i);
	}
}

