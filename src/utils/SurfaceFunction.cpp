#include "SurfaceFunction.hpp"

#include <string>
#include <stdexcept>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

inline void swap(SurfaceFunction & left, SurfaceFunction & right)
{
    using std::swap;
    swap(left.name, right.name);
    swap(left.nPoints, right.nPoints); // This is maybe redundant as nPoints must be the same...
    swap(left.values, right.values);
}

inline void SurfaceFunction::swap(SurfaceFunction & other)
{
    using std::swap;
    swap(this->name, other.name);
    swap(this->nPoints, other.nPoints); // This is maybe redundant
    swap(this->values, other.values);
}

SurfaceFunction & SurfaceFunction::operator=(SurfaceFunction other)
{
    if (this == &other) // Check for self-assignment
        throw std::runtime_error("Are you trying to self-assign?");
    // No check on dimensions! This is assignment not comparison!!

    ::swap(*this, other);
    return *this;
}

double SurfaceFunction::operator*(const SurfaceFunction & other)
{
	if (this->nPoints != other.nPoints)
		throw std::runtime_error("Incoherent dimensions of left and right operands!");
	return this->values.dot(other.values);
}

SurfaceFunction & SurfaceFunction::operator+=(const SurfaceFunction & other)
{
	if (this->nPoints != other.nPoints)
		throw std::runtime_error("Incoherent dimensions of left and right operands!");
	this->name += "+" + other.name;
	this->values = this->values + other.values;
    	return *this;
}

SurfaceFunction & SurfaceFunction::operator-=(const SurfaceFunction & other)
{
	if (this->nPoints != other.nPoints)
		throw std::runtime_error("Incoherent dimensions of left and right operands!");
	this->name += "-" + other.name;
	this->values = this->values - other.values;
    	return *this;
}

SurfaceFunction & SurfaceFunction::operator*=(double scaling)
{
#if defined (HAS_CXX11)
	// The C++11 way
	// As from the standard this would produce the same result at std::sprintf(buf, "%f", value)
	// meaning that 2.5 will be represented as 2.500000
	std::string tmp = std::to_string(scaling);
	tmp.erase(tmp.find_last_not_of('0') + 1, std::string::npos);
	this->name = tmp + "*" + this->name; 
#else
	std::ostringstream sstream;
	sstream << scaling;
	std::string scalingAsString = sstream.str();
	this->name = scalingAsString + "*" + this->name;
#endif
	this->values *= scaling;
        return *this;
}

void SurfaceFunction::setValues(double * values_) 
{
	if (!allocated)
		throw std::runtime_error("Surface function not allocated!");
	// Zero out any previous value
	values.setZero();

	for (int i = 0; i < nPoints; ++i) 
	{
		values(i) = values_[i];
		}
}

void SurfaceFunction::getValues(double * values_) 
{
	for (int i = 0; i < nPoints; ++i) 
	{
		values_[i] = values(i);
	}
}

void SurfaceFunction::clear() 
{
	values.setZero();
}

std::ostream & SurfaceFunction::printObject(std::ostream & os) 
{
	os << "Surface Function " << name << std::endl;
	if (!allocated) 
		throw std::runtime_error("Surface function not allocated!");
	
	os << values.transpose();
	
	return os;
}
