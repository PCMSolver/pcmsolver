#include "Sphere.hpp"

#include <ostream>
#include <stdexcept>

#include "Config.hpp"

std::ostream & Sphere::printObject(std::ostream & os) 
{
	os << "Sphere radius " << radius_ << std::endl;
	os << "Sphere center\n" << center_;

	return os;
}

inline void swap(Sphere & left, Sphere & right)
{
    using std::swap;
    swap(left.center_, right.center_);
    swap(left.radius_, right.radius_);
    swap(left.colour_, right.colour_);
}

inline void Sphere::swap(Sphere & other)
{
    using std::swap;
    swap(this->center_, other.center_);
    swap(this->radius_, other.radius_);
    swap(this->colour_, other.colour_);
}

Sphere & Sphere::operator=(Sphere other)
{
    if (this == &other) // Check for self-assignment
        throw std::runtime_error("Are you trying to self-assign?");

    ::swap(*this, other);
    return *this;
}
