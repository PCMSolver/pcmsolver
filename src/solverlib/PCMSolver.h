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

#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "GreensFunctionSum.h"
#include "MetalSphere.h"

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
    virtual void buildPCMMatrix();
    virtual bool readCavity(string &filename);
 private:
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
