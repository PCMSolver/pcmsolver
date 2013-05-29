#include <fstream>
#include <stdexcept>

#include "SurfaceFunction.h"

SurfaceFunction & SurfaceFunction::operator=(const SurfaceFunction & other)
{
    if (this == &other) // Check for self-assignment
        throw std::runtime_error("Are you trying to self-assign?");
    if (this->nPoints != other.nPoints) // Check if dimensions match
        throw std::runtime_error("Incoherent dimensions!");

    // Now swap this and other using the Copy-constructor and a non-throwing swap
    SurfaceFunction(other).swap(*this);
    // Old resources are released with the destruction of the temporary above
    return *this;
}

SurfaceFunction operator+(const SurfaceFunction & left, const SurfaceFunction & right)
{
    if (left.nPoints != right.nPoints)
        throw std::runtime_error("Incoherent dimensions of left and right operands!");
    std::string resultName = left.name + "+" + right.name;
    SurfaceFunction result(resultName, left.nPoints);
    result.values = left.values + right.values;
    return result;
}

SurfaceFunction operator-(const SurfaceFunction & left, const SurfaceFunction & right)
{
    if (left.nPoints != right.nPoints)
        throw std::runtime_error("Incoherent dimensions of left and right operands!");
    std::string resultName = left.name + "-" + right.name;
    SurfaceFunction result(resultName, left.nPoints);
    result.values = left.values - right.values;
    return result;
}

SurfaceFunction & SurfaceFunction::operator+=(const SurfaceFunction & other)
{
    *this = *this + other;
    return *this;
}

SurfaceFunction & SurfaceFunction::operator-=(const SurfaceFunction & other)
{
    *this = *this - other;
    return *this;
}

inline void SurfaceFunction::swap(SurfaceFunction & other)
{
    if (this->nPoints != other.nPoints) // Check if dimensions match
        throw std::runtime_error("Incoherent dimensions!");
    using std::swap;
    swap(this->name, other.name);
    swap(this->allocate(nPoints, other.nPoints); // This is maybe redundant
    swap(this->values, other.values);
}

inline void SurfaceFunction::swap(SurfaceFunction & left, SurfaceFunction & right)
{
    if (left.nPoints != right.nPoints) // Check if dimensions match
        throw std::runtime_error("Incoherent dimensions!");
    using std::swap;
    swap(left.name, right.name);
    swap(left.nPoints, right.nPoints); // This is maybe redundant
    swap(left.values, right.values);
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

std::ostream & operator<<(std::ostream & os, SurfaceFunction & sf) 
{
	return sf.printObject(os);
}

std::ostream & SurfaceFunction::printObject(std::ostream & os) {
	os << "Surface Function " << name << std::endl;
	if (!allocated) 
		throw std::runtime_error("Surface function not allocated!");
	
	os << values.transpose() << std::endl;
	
	return os;
}

