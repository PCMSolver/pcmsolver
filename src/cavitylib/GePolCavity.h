#ifndef GEPOLCAVITY
#define GEPOLCAVITY

#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "Getkw.h"
#include "Cavity.h"
#include "Atom.h"
#include "Sphere.h"

//class Getkw;

/*

  C++ inteface and wrapper for GePol cavity class
  written by Krzysztof Mozgawa, 2011

*/

class Getkw;

class GePolCavity : public Cavity 
{
 public:
    GePolCavity(){}
    GePolCavity(const Section & cavity);
    ~GePolCavity(){};
    void makeCavity(int maxts, int lwork);
    void makeCavity();
    void writeOutput(string &filename);
    VectorXd & getTessRadius(){return tessRadius;}
    int getNSpheres(){return nSpheres;}
    void setNSpheres(int n){nSpheres = n;}
    Matrix<double, 3, Dynamic> & getTessSphereCenter(){return tessSphereCenter;}
    double getTessRadius(int i){return tessRadius(i);}
    void setMaxAddedSpheres(bool add = true, int maxAdd = 100);
    double getProbeRadius() { return probeRadius; };
    void setProbeRadius( double probeRadius );
    vector<Sphere> & getSpheres(){ return spheres; }
    int getMode(){return mode;}
    enum SphereMode {Explicit, Atoms, Implicit};

 protected:
    virtual ostream & printObject(ostream & os);
        
 private:
    SphereMode mode;
    void setMode(const std::string & mode);
    void setMode(int mode);
    int nSpheres;
    int maxAddedSpheres;
    int addedSpheres;
    bool addSpheres;
    double probeRadius;
    Matrix<double, 3, Dynamic> tessSphereCenter;
    VectorXd tessRadius;
    vector<Sphere> spheres;
};


#endif
