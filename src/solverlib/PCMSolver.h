/*! \file PCMSolver.h 
\brief PCM solver
*/


#ifndef PCMSOLVER
#define PCMSOLVER

#include <string>
#include <vector>
#include <iostream>
#include <complex>

using namespace std;

class PCMSolver{
 public:
    PCMSolver(GreensFunction &gfi, GreensFunction &gfo){
        greenInside = &gfi;
        greenOutside = &gfo;
    };
    ~PCMSolver(){};
    GreensFunction &getGreenInside();
    GreensFunction &getGreenOutside();    
    int getCavitySize() const {return cavitySize;};
    const MatrixXd& getPCMMatrix() const {return PCMMatrix;};
    virtual void buildPCMMatrix(GePolCavity cav);
    virtual bool readCavity(string &filename);
    virtual double compDiagonalElementSoper(GreensFunction *green, int i, GePolCavity cav);
    virtual double compDiagonalElementDoper(GreensFunction *green, int i, GePolCavity cav);
 private:
    static const double factor = 1.0694;
    int cavitySize;
    GreensFunction *greenInside;
    GreensFunction *greenOutside;
    VectorXd areaTess;
    VectorXd radiusTess;
    Matrix<double, Dynamic, 3> centerSphereTess;
    Matrix<double, Dynamic, 3> centerTess;
    Matrix<double, Dynamic, 3> normalTess;
    MatrixXd PCMMatrix;
};
#endif
