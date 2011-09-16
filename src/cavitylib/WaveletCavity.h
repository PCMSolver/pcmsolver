#ifndef WAVELETCAVITY
#define WAVELETCAVITY

#include <iostream>
#include <string>

#include <Eigen/Dense>

#include "Getkw.h"
#include "Cavity.h"

class Getkw;

class WaveletCavity : public Cavity {
 public:
    WaveletCavity(){}
    //    WaveletCavity(string &filename);
    WaveletCavity(const Getkw & input, const string path = "Cavity");
    WaveletCavity(const Section & cavity);
    ~WaveletCavity(){};
    void makeCavity();
    VectorXd & getTessRadius(){return tessRadius;};
    VectorXd & getSphereRadius(){return sphereRadius;};
    int getNSpheres(){return nSpheres;};
    Matrix<double, 3, Dynamic> & getSphereCenter(){return sphereCenter;};
    Matrix<double, 3, Dynamic> & getTessSphereCenter(){return tessSphereCenter;};
    double getTessRadius(int i){return tessRadius(i);};

    friend std::ostream& operator<<(std::ostream &o, const WaveletCavity &c);

 private:
    void writeInput(string &fileName);
    int nSpheres;
    Matrix<double, 3, Dynamic> sphereCenter;
    VectorXd sphereRadius;
    Matrix<double, 3, Dynamic> tessSphereCenter;
    VectorXd tessRadius;
    int patchLevel;
    double probeRadius;
    double coarsity;
};

#endif
