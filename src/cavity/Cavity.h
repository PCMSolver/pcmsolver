#ifndef CAVITY
#define CAVITY

#include <iostream>
#include <string>
#include <map>

#include "Config.h"

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;


/** 

A basic cavity class
written by Krzysztof Mozgawa, 2011


*/

class SurfaceFunction;

typedef std::pair< std::string, SurfaceFunction * > SurfaceFunctionPair;
typedef std::map< std::string, SurfaceFunction * > SurfaceFunctionMap;

class Cavity
{
 public:
    Cavity(){built = false; nElements = 0;}
    ~Cavity(){}
    virtual void makeCavity() = 0;
    virtual void writeOutput(string &filename);
    virtual Matrix3Xd & getElementCenter(){return elementCenter;}
    virtual Vector3d getElementCenter(int i){return elementCenter.col(i);}
    virtual Matrix3Xd & getElementNormal(){return elementNormal;}
    virtual Vector3d getElementNormal(int i){return elementNormal.col(i);}
    virtual VectorXd & getElementArea(){return elementArea;}
    virtual double getElementArea(int i){return elementArea(i);}
    virtual int size(){return nElements;}
    bool isBuilt(){return built;}
    double compPolarizationEnergy();
    double compPolarizationEnergy(const std::string & potential, 
                                  const std::string & charge);
    void appendNewFunction(const std::string & name);
    void setFunction(const std::string & name, double * values);
    SurfaceFunction & getFunction(const std::string & name);
    bool functionExists(const std::string & name) { 
        return (functions.count(name) == 1);
    }
    enum chargeType{Nuclear, Electronic};

    friend std::ostream& operator<<(std::ostream & o, Cavity & c);
    
 protected:
    virtual ostream & printObject(ostream & os);
    int nElements;
    bool built;
    Matrix3Xd elementCenter;
    Matrix3Xd elementNormal;
    VectorXd elementArea;
    SurfaceFunctionMap functions;
};


#endif
