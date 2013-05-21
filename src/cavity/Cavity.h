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
    Cavity(){built = false; nTess = 0;}
    ~Cavity(){}
    virtual void makeCavity() = 0;
    virtual void writeOutput(string &filename);
    virtual Matrix<double, 3, Dynamic> & getTessCenter(){return tessCenter;}
    virtual Vector3d getTessCenter(int i){return tessCenter.col(i);}
    virtual Matrix<double, 3, Dynamic> & getTessNormal(){return tessNormal;}
    virtual Vector3d getTessNormal(int i){return tessNormal.col(i);}
    virtual VectorXd & getTessArea(){return tessArea;}
    virtual double getTessArea(int i){return tessArea(i);}
    virtual int size(){return nTess;}
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
    int nTess;
    bool built;
    Matrix<double, 3, Dynamic> tessCenter;
    Matrix<double, 3, Dynamic> tessNormal;
    VectorXd tessArea;
    SurfaceFunctionMap functions;
};


#endif
