/*! \file GreensFunctionInterface.h 
\brief Interface for common operations for Green´s functions.
*/

#ifndef GREENSFUNCTIONINTERFACE
#define GREENSFUNCTIONINTERFACE

/** Green´s function interface

A generic green´s function to reprensent the electrostatic potential for a given environment

*/

class Section;

class GreensFunctionInterface
{
 public:
    virtual ~GreensFunctionInterface() {};
    virtual double evalf(Vector3d &p1, Vector3d &p2) = 0;
    virtual double evald(Vector3d &direction, Vector3d &p1, Vector3d &p2) = 0;
    virtual double derivativeSource(Vector3d &direction, Vector3d &p1, Vector3d &p2) = 0;
    virtual double derivativeProbe(Vector3d &direction, Vector3d &p1, Vector3d &p2) = 0;
    virtual Vector3d gradientSource(Vector3d &p1, Vector3d &p2) = 0;
    virtual Vector3d gradientProbe(Vector3d &p1, Vector3d &p2) = 0;
    virtual void gradientSource(Vector3d &gradient, Vector3d &p1, Vector3d &p2) = 0;
    virtual void gradientProbe(Vector3d &gradient, Vector3d &p1, Vector3d &p2) = 0;
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
