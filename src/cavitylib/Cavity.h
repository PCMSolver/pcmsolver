#ifndef CAVITY
#define CAVITY

#include <iostream>
#include <string>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;


/** 

A basic cavity class
written by Krzysztof Mozgawa, 2011


*/

class Cavity
{
 public:
    Cavity(){}
    ~Cavity(){}
    virtual void makeCavity(int, int) = 0; // not nice, needs fix
    virtual void writeOutput(string &filename);
    virtual Matrix<double, Dynamic, 3> & getTessCenter(){return tessCenter;}
    virtual Matrix<double, Dynamic, 3> & getTessNormal(){return tessNormal;}
    virtual VectorXd & getTessArea(){return tessArea;}
    virtual int size(){return nTess;}
 protected:
    int nTess;
    Matrix<double, Dynamic, 3> tessCenter;
    Matrix<double, Dynamic, 3> tessNormal;
    VectorXd tessArea;
    double averageArea;
};


#endif
