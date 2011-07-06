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
    virtual void makeCavity() = 0;
    virtual void writeOutput(string &filename);
    virtual int getNTess(){ return nTess;}
    virtual Matrix<double, Dynamic, 3> & getTessCenter(){return tessCenter;}
    virtual Matrix<double, Dynamic, 3> & getTessNormal(){return tessNormal;}
    virtual VectorXd & getTessArea(){return tessArea;}
    virtual double getTessArea(int i){return tessArea(i);}
    virtual double getTessCenter(int i, int j){return tessCenter(i,j);}
    virtual int size(){return nTess;}
 protected:
    int nTess;
    Matrix<double, Dynamic, 3> tessCenter;
    Matrix<double, Dynamic, 3> tessNormal;
    VectorXd tessArea;
    double averageArea;
};


#endif
