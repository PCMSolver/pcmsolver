#ifndef GEPOLCAVITY
#define GEPOLCAVITY

#include <iostream>
#include <string>
#include <vector>

#include <Config.h>

#include "CavityOfSpheres.h"


/*

  C++ inteface and wrapper for GePol cavity class
  written by Krzysztof Mozgawa, 2011

*/

class GePolCavity : public CavityOfSpheres 
{
 public:
    GePolCavity(){}
    GePolCavity(double _area, const std::vector<Sphere> & _spheres, bool _addSpheres = false, double _probeRadius = 0.0) : 
     averageArea(_area), addSpheres(_addSpheres), probeRadius(_probeRadius) 
           {
		// Initialize data members inherited from CavityOfSpheres...
		spheres = _spheres;
                nSpheres = spheres.size();
		sphereCenter.resize(Eigen::NoChange, nSpheres);
		sphereRadius.resize(nSpheres);
		for (int i = 0; i < nSpheres; ++i) {
			sphereCenter.col(i) = spheres[i].getSphereCenter();
			sphereRadius(i) = spheres[i].getSphereRadius();
		}
		// ...and build the cavity!
		makeCavity(10000, 10000000);
           }
    ~GePolCavity(){}
    void makeCavity(int maxts, int lwork);
    void makeCavity();
    void writeOutput(string &filename);
    void setMaxAddedSpheres(bool add = true, int maxAdd = 100);
    double getProbeRadius() { return probeRadius; }
    void setProbeRadius( double probeRadius );

 protected:
    virtual ostream & printObject(ostream & os);
        
 private:
    double averageArea;
    bool addSpheres;
    double probeRadius;
    int maxAddedSpheres;
    int addedSpheres;
};


#endif
