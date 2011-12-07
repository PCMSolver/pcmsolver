/*! \file GreensFunction.h 
\brief Abstract base class for the Green´s function generator.
*/

#ifndef GREENSFUNCTION
#define GREENSFUNCTION

/** Green´s function Abstract base class 

A generic green´s function to reprensent the electrostatic potential for a given environment

*/
class GreensFunction
{
 public:
    GreensFunction(){delta = 1.0e-4;}
    virtual ~GreensFunction(){};
    virtual double evalf(Vector3d &p1, Vector3d &p2) = 0;
    virtual double evald(Vector3d &direction, Vector3d &p1, Vector3d &p2) = 0;
    virtual double derivative(Vector3d &direction, Vector3d &p1, Vector3d &p2);
    virtual void gradient(Vector3d &gradient, Vector3d &p1, Vector3d &p2);
    void setDelta(double value);
    double getDelta(){return delta;}
    bool isUniform(){ return uniformFlag; }
    GreensFunction * allocateGreensFunction(const Section &green);
    GreensFunction * allocateGreensFunction(double dielConst);
    GreensFunction * allocateGreensFunction();
 protected:
    double delta;
	bool uniformFlag;	
};
#endif
