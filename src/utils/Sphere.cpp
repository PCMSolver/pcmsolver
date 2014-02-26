#include "Sphere.hpp"

#include <ostream>
#include <stdexcept>

#include "Config.hpp"

std::ostream & Sphere::printObject(std::ostream & os) 
{
	os << "Sphere radius " << sphereRadius_ << std::endl;
	os << "Sphere center\n" << sphereCenter_;

	return os;
}

inline void swap(Sphere & left, Sphere & right)
{
    using std::swap;
    swap(left.sphereCenter_, right.sphereCenter_);
    swap(left.sphereRadius_, right.sphereRadius_);
    swap(left.sphereColour_, right.sphereColour_);
}

inline void Sphere::swap(Sphere & other)
{
    using std::swap;
    swap(this->sphereCenter_, other.sphereCenter_);
    swap(this->sphereRadius_, other.sphereRadius_);
    swap(this->sphereColour_, other.sphereColour_);
}

Sphere & Sphere::operator=(Sphere other)
{
    if (this == &other) // Check for self-assignment
        throw std::runtime_error("Are you trying to self-assign?");

    ::swap(*this, other);
    return *this;
}
