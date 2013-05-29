#ifndef SURFACE_FUNCTION
#define SURFACE_FUNCTION

#include <iostream>
#include <string>

#include "Config.h"

#include <Eigen/Dense>

/** 

A basic surface function class
written by L. Frediani 2012

*/

class SurfaceFunction
{
 public:
    SurfaceFunction(const std::string & name_) : name(name_)
	{
		allocated = false;
	}
    SurfaceFunction(const std::string & name_, int nPoints_) : name(name_), nPoints(nPoints_) 
	{
		values.resize(nPoints);
		allocated = true;
	}							       
    SurfaceFunction(const std::string & name_, int nPoints_, double * values_) : name(name_), nPoints(nPoints_)
	{
		values.resize(nPoints);
		allocated = true;
		for (int i = 0; i < nPoints; ++i)
		{
			values(i) = values_[i];
		}
	}
    ~SurfaceFunction()
    {
        allocated = false;
    }

    /// Copy constructor
    SurfaceFunction(const SurfaceFunction & other) : name(other.name), nPoints(other.nPoints), values(other.values)
    {
        allocated = true;
    }

    /// Assignment operator
    SurfaceFunction & operator=(const SurfaceFunction & other);
    /// Addition operator
    friend SurfaceFunction operator+(const SurfaceFunction & left, const SurfaceFunction & right);
    /// Subtraction operator
    friend SurfaceFunction operator-(const SurfaceFunction & left, const SurfaceFunction & right);
    /// Addition-assignment operator
    SurfaceFunction & operator+=(const SurfaceFunction & other);
    /// Subtraction-assignment operator
    SurfaceFunction & operator-=(const SurfaceFunction & other);

    inline void swap(SurfaceFunction & left, SurfaceFunction & right);
    inline void swap(SurfaceFunction & other);
    void setValue(int index_, double value_) { values(index_) = value_; }
    double getValue(int index_) {return values(index_);}
    Eigen::VectorXd & getVector(){ return values; }
    void allocate(int nPoints_){ values.resize(nPoints_); }
    bool isAllocated() { return allocated; }
    void clear();

    void setValues(double * value);
    void getValues(double * value);

    friend std::ostream & operator<<(std::ostream & o, SurfaceFunction & s);

 private:
    virtual std::ostream & printObject(std::ostream & os);
    std::string name;
    int nPoints;
    Eigen::VectorXd values;
    bool allocated;
};

#endif
