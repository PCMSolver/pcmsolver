#ifndef GEPOLCAVITY
#define GEPOLCAVITY

#include <iostream>
#include <string>

#include <Eigen/Dense>

/*

C++ inteface and wrapper for GePol cavity class
written by Krzysztof Mozgawa, 2011

*/
class GePolCavity : public Cavity {
 public:
    GePolCavity(string &filename){
        readInput(filename);
    };
    ~GePolCavity(){};
    void makeCavity();
    void writeOutput(string &filename);
    VectorXd & getTessRadius(){return tessRadius;};
    VectorXd & getSphereRadius(){return sphereRadius;};
    int getNSpheres(){return nSpheres;};
    Matrix<double, Dynamic, 3> & getSphereCenter(){return sphereCenter;};
    Matrix<double, Dynamic, 3> & getTessSphereCenter(){return tessSphereCenter;};
    double getTessRadius(int i){return tessRadius(i);};
 private:
    bool readInput(string &filename);
    int nSpheres;
    VectorXd sphereRadius;
    Matrix<double, Dynamic, 3> sphereCenter;
    Matrix<double, Dynamic, 3> tessSphereCenter;
    VectorXd tessRadius;
    // Variables needed for communication with pedra cavity
    // if one would require more thatn 300 spheres it needs to be changed here
    double nesfp;
    double xe[300], ye[300], ze[300], rin[300];

};


#endif
