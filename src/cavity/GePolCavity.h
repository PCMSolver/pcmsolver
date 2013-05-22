#ifndef GEPOLCAVITY
#define GEPOLCAVITY

#include <iostream>
#include <string>
#include <vector>

#include <Config.h>

#include <Eigen/Dense>

#include "Getkw.h"
#include "CavityOfSpheres.h"


/*

  C++ inteface and wrapper for GePol cavity class
  written by Krzysztof Mozgawa, 2011

*/

class GePolCavity : public CavityOfSpheres 
{
 public:
    GePolCavity(){}
    GePolCavity(const Section & cavity);
    GePolCavity(const Section & cavity, const vector<Sphere> & _spheres);
    GePolCavity(double _area, const std::vector<Sphere> & _spheres, bool _addSpheres = false, double _probeRadius = 0.0) : 
     averageArea(_area), addSpheres(_addSpheres), probeRadius(_probeRadius) 
           {
		spheres = _spheres;
                nSpheres = spheres.size();
           }

    ~GePolCavity(){};
    void makeCavity(int maxts, int lwork);
    void makeCavity();
    void writeOutput(string &filename);
    void setMaxAddedSpheres(bool add = true, int maxAdd = 100);
    double getProbeRadius() { return probeRadius; };
    void setProbeRadius( double probeRadius );
    /*
    VectorXd & getElementRadius(){return elementRadius;}
    int getNSpheres(){return nSpheres;}
    void setNSpheres(int n){nSpheres = n;}
    Matrix3Xd & getElementSphereCenter(){return elementSphereCenter;}
    double getElementRadius(int i){return elementRadius(i);}
    vector<Sphere> & getSpheres(){ return spheres; }
    int getMode(){return mode;}
    enum SphereMode {Explicit, Atoms, Implicit};*/

 protected:
    virtual ostream & printObject(ostream & os);
        
 private:
    double averageArea;
    //vector<Sphere> spheres;
    bool addSpheres;
    double probeRadius;
    //SphereMode mode;
    //void setMode(const std::string & mode);
    //void setMode(int mode);
    //int nSpheres;
    int maxAddedSpheres;
    int addedSpheres;
    //Matrix3Xd elementSphereCenter;
    //VectorXd elementRadius;
};


#endif
