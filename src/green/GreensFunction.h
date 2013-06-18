/*! \file GreensFunction.h 
\brief Abstract base class for the Green´s function generator.
*/

#ifndef GREENSFUNCTION
#define GREENSFUNCTION

/** Green´s function Abstract base class 

A generic green´s function to reprensent the electrostatic potential for a given environment

*/

#include "GreensFunctionInterface.h"

class Section;

template<class T>
class GreensFunction: public GreensFunctionInterface
{
 public:
    GreensFunction(){delta = 1.0e-4;}
    virtual ~GreensFunction(){};

    // From GreensFunctionInterface
    virtual double evalf(Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    virtual double evald(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2) = 0;
    virtual double derivativeSource(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    virtual double derivativeProbe(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    virtual Eigen::Vector3d gradientSource(Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    virtual Eigen::Vector3d gradientProbe(Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    virtual double getDielectricConstant();
    virtual void gradientSource(Eigen::Vector3d &gradient, Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    virtual void gradientProbe(Eigen::Vector3d &gradient, Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    void setDelta(double value);
    double getDelta(){return delta;}
    GreensFunction<T> * allocateGreensFunction(const Section &green);
    GreensFunction<T> * allocateGreensFunction(double dielConst);
    GreensFunction<T> * allocateGreensFunction();
    virtual double compDiagonalElementS(double area) = 0;
    virtual double compDiagonalElementD(double area, double radius) = 0;
 protected:
    virtual T evalGreensFunction(T * source, T * probe) = 0;
    std::ostream & printObject(std::ostream & os);
    double delta;
    bool uniformFlag;
};


#endif
