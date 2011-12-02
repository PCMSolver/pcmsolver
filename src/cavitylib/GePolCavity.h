#ifndef GEPOLCAVITY
#define GEPOLCAVITY

#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "Getkw.h"
#include "Cavity.h"
#include "Atom.h"
//class Getkw;

/*

C++ inteface and wrapper for GePol cavity class
written by Krzysztof Mozgawa, 2011

*/

class Getkw;

class GePolCavity : public Cavity {
 public:
    GePolCavity(){}
    GePolCavity(string &filename){
        readInput(filename);
    }
    GePolCavity(const Getkw & Input, const string path = "Cavity");
    GePolCavity(const Section & cavity);
    ~GePolCavity(){};
    void makeCavity(int maxts, int lwork);
    void makeCavity();
    void writeOutput(string &filename);
    VectorXd & getTessRadius(){return tessRadius;}
    VectorXd & getSphereRadius(){return sphereRadius;}
    int getNSpheres(){return nSpheres;}
    void setNSpheres(int n){nSpheres = n;}
    VectorXd & getCharges(){return charges;}
    Matrix<double, 3, Dynamic> & getSphereCenter(){return sphereCenter;}
    Matrix<double, 3, Dynamic> & getTessSphereCenter(){return tessSphereCenter;}
    double getTessRadius(int i){return tessRadius(i);}


    friend std::ostream& operator<<(std::ostream &o, const GePolCavity &c);

 private:
    bool readInput(string &filename);
    int nSpheres;
    VectorXd sphereRadius;
    VectorXd charges;
    Matrix<double, 3, Dynamic> sphereCenter;
    Matrix<double, 3, Dynamic> tessSphereCenter;
    VectorXd tessRadius;
   
    // Variables needed for communication with pedra cavity
    // if one would require more thatn 300 spheres it needs to be changed here
};


#endif
