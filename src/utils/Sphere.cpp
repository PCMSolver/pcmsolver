#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "Sphere.h"

/*

  Methods for Sphere class
  written by Roberto Di Remigio, 2011

*/

ostream & operator<<(ostream & os, Sphere & sphere) {
	return sphere.printObject(os);
} 

ostream & Sphere::printObject(ostream & os) {
	os << "Sphere radius " << sphereRadius << endl;
	os << "Sphere center\n" << sphereCenter;
	return os;
}

inline void swap(Sphere & left, Sphere & right)
{
    using std::swap;
    swap(left.sphereCenter, right.sphereCenter);
    swap(left.sphereRadius, right.sphereRadius);
    swap(left.sphereColour, right.sphereColour);
}

inline void Sphere::swap(Sphere & other)
{
    using std::swap;
    swap(this->sphereCenter, other.sphereCenter);
    swap(this->sphereRadius, other.sphereRadius);
    swap(this->sphereColour, other.sphereColour);
}

Sphere & Sphere::operator=(Sphere other)
{
    if (this == &other) // Check for self-assignment
        throw std::runtime_error("Are you trying to self-assign?");

    ::swap(*this, other);
    return *this;
}
