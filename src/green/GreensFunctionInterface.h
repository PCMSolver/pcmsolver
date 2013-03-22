/*! \file GreensFunctionInterface.h 
\brief Interface for common operations for Green´s functions.
*/

#ifndef GREENSFUNCTIONINTERFACE
#define GREENSFUNCTIONINTERFACE

#include <Eigen/Dense>

/** Green´s function interface

A generic green´s function to reprensent the electrostatic potential for a given environment

*/

class Section;

class GreensFunctionInterface
{
 public:
    virtual ~GreensFunctionInterface() {};
    virtual double evalf(Eigen::Vector3d &p1, Eigen::Vector3d &p2) = 0;
    virtual double evald(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2) = 0;
    virtual double derivativeSource(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2) = 0;
    virtual double derivativeProbe(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2) = 0;
    virtual Eigen::Vector3d gradientSource(Eigen::Vector3d &p1, Eigen::Vector3d &p2) = 0;
    virtual Eigen::Vector3d gradientProbe(Eigen::Vector3d &p1, Eigen::Vector3d &p2) = 0;
    virtual void gradientSource(Eigen::Vector3d &gradient, Eigen::Vector3d &p1, Eigen::Vector3d &p2) = 0;
    virtual void gradientProbe(Eigen::Vector3d &gradient, Eigen::Vector3d &p1, Eigen::Vector3d &p2) = 0;
    virtual double getDielectricConstant() = 0;
    virtual double compDiagonalElementS(double area) = 0;
    virtual double compDiagonalElementD(double area, double radius) = 0;
    GreensFunctionInterface * allocateGreensFunctionInterface(const Section &green);
    GreensFunctionInterface * allocateGreensFunctionInterface(double dielConst, const std::string greenDer = "Derivative");
    GreensFunctionInterface * allocateGreensFunctionInterface(const std::string greenDer = "Derivative");
    bool isUniform(){ return uniformFlag; }
    enum derivativeType{Numerical, Directional, Gradient, Hessian};
    friend std::ostream & operator<<(std::ostream & os, GreensFunctionInterface &gf){
        return gf.printObject(os);
    };

 protected:
    virtual std::ostream & printObject(std::ostream & os) = 0;
    bool uniformFlag;
};


#endif
