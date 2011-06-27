#include <iostream>
#include <ostream>
#include <cmath>
#include <cstdlib>

using namespace std;

#include "GreensFunction.h"
#include "PlanarInterface.h"

PlanarInterface::PlanarInterface(double eps1, double eps2, double pos, double width)
{
    this->eps1 = eps1;
    this->eps2 = eps2;
    this->pos = pos;
    this->width = width;
    computed = false;
}

double PlanarInterface::evalf(double* p1, double* p2) {
    if (computed) {
	cout << "method not implemented yet" << endl;
	exit(-1);
    }
    else {
	cout << "compute green´s function before using it!" << endl;
	exit(-1);
    }
}
