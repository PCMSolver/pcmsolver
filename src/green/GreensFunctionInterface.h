/*! \file GreensFunctionInterface.h 
\brief Interface for common operations for Green´s functions.
*/

#ifndef GREENSFUNCTIONINTERFACE
#define GREENSFUNCTIONINTERFACE

/** Green´s function interface

A generic green´s function to reprensent the electrostatic potential for a given environment

*/

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
    virtual void compDiagonalElementS(double area) = 0;
    virtual void compDiagonalElementD(double area, double radius) = 0;
    GreensFunctionInterface * allocateGreensFunctionInterface(const Section &green);
    GreensFunctionInterface * allocateGreensFunctionInterface(double dielConst);
    GreensFunctionInterface * allocateGreensFunctionInterface();

    bool isUniform(){ return uniformFlag; }

    enum derivativeType{Numerical, Directional, Gradient, Hessian};
    enum greensFunctionType{Vacuum, UniformDielectric, MetalSphere};

 protected:
    bool uniformFlag;
};


#endif
