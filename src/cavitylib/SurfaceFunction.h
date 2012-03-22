#ifndef SURFACE_FUNCTION
#define SURFACE_FUNCTION

#include <iostream>
#include <string>
#include <Eigen/Dense>

using namespace Eigen;

/** 

A basic surface function class
written by L. Frediani 2012

*/

class SurfaceFunction
{
 public:
    SurfaceFunction(std::string & name);
    SurfaceFunction(std::string & name, int nPoints);
    SurfaceFunction(std::string & name, int nPoints, double * values);
    ~SurfaceFunction(){}

    void setValue(int index, double value) {values(index) = value;}
    double getValue(int index) {return values(index);}
    VectorXd getVector(){return values;}
    void allocate(int nPoints){values.resize(nPoints);}
    bool isAllocated(){return allocated;}

    void setValues(double * value);
    void getValues(double * value);
 private:
    bool allocated;
    std::string name;
    VectorXd values;
};

#endif
