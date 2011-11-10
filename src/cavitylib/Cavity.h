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
    Cavity(){isBuilt = false; nTess = 0;}
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

    virtual void initPotChg();

    VectorXd & getChg(int type);
    VectorXd & getPot(int type);
    double getChg(int type, int i);
    double getPot(int type, int i);

    double compPolarizationEnergy();

    enum chargeType{Nuclear, Electronic};

    friend std::ostream& operator<<(std::ostream &o, const Cavity &c);

 protected:
    int nTess;
    bool isBuilt;
    Matrix<double, 3, Dynamic> tessCenter;
    Matrix<double, 3, Dynamic> tessNormal;
    VectorXd tessArea;
    double averageArea;
    VectorXd nuclearPotential;
    VectorXd nuclearCharge;
    VectorXd electronicPotential;
    VectorXd electronicCharge;
};


#endif
