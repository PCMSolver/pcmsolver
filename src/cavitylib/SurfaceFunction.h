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
    SurfaceFunction(const std::string & name);
    SurfaceFunction(const std::string & name, int nPoints);
    SurfaceFunction(const std::string & name, int nPoints, double * values);
    ~SurfaceFunction(){}

    void setValue(int index, double value) {values(index) = value;}
    double getValue(int index) {return values(index);}
    VectorXd & getVector(){return values;}
    void allocate(int nPoints){values.resize(nPoints);}
    bool isAllocated(){return allocated;}
    void clear();

    void setValues(double * value);
    void getValues(double * value);

    friend std::ostream & operator<<(std::ostream & o, SurfaceFunction & s);

 private:
    virtual std::ostream & printObject(std::ostream & os);
    bool allocated;
    std::string name;
    VectorXd values;
};

#endif
