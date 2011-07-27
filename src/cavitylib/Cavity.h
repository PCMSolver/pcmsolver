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
    virtual Vector3d getTessCenter(int i){return tessCenter.row(i);}

    virtual Matrix<double, Dynamic, 3> & getTessNormal(){return tessNormal;}
    virtual Vector3d getTessNormal(int i){return tessNormal.row(i);}

    virtual VectorXd & getTessArea(){return tessArea;}
    virtual double getTessArea(int i){return tessArea(i);}

    virtual int size(){return nTess;}

    friend std::ostream& operator<<(std::ostream &o, const Cavity &c);
 protected:
    int nTess;
    Matrix<double, Dynamic, 3> tessCenter;
    Matrix<double, Dynamic, 3> tessNormal;
    VectorXd tessArea;
    double averageArea;
};


#endif
