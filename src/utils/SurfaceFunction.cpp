#include "SurfaceFunction.hpp"

#include <map>
#include <string>
#include <stdexcept>

#include "Config.hpp"

#include <Eigen/Dense>

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
	this->name = std::to_string(scaling) + "*" + this->name;
	this->values *= scaling;
        return *this;
}

void SurfaceFunction::setValues(double * values_) 
{
	if (!allocated)
		throw std::runtime_error("Surface function not allocated!");

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

bool SurfaceFunction::Register()
{
	registered = SurfaceFunction::TheMap().insert(SurfaceFunctionMap::value_type(name, this)).second;
	return registered;
}

bool SurfaceFunction::unRegister()
{
	registered = SurfaceFunction::TheMap().erase(name);
	return registered;
}

std::ostream & SurfaceFunction::printObject(std::ostream & os) {
	os << "Surface Function " << name << std::endl;
	if (!allocated) 
		throw std::runtime_error("Surface function not allocated!");
	
	os << values.transpose();
	
	return os;
}
