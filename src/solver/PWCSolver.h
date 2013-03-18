/*! \file PWCSolver.h 
\brief PCM solver
*/


#ifndef PWCSOLVER_H_
#define PWCSOLVER_H_

#include <string>
#include <vector>
#include <iostream>
#include <complex>

extern "C"{
#include "vector3.h"
#include "sparse2.h"
#include "intvector.h"
#include "basis.h"
}

using namespace std;

class PWCSolver : public WEMSolver {
 public:
    PWCSolver(GreensFunctionInterface & gfi, GreensFunctionInterface & gfo);
    PWCSolver(GreensFunctionInterface * gfi, GreensFunctionInterface * gfo);
    PWCSolver(const Section & solver);
    ~PWCSolver();
 private: 
    virtual void initInterpolation();
    virtual void constructWavelets();
    virtual void constructSi();
    virtual void constructSe();
    virtual void solveFirstKind(const VectorXd & potential, VectorXd & charge);
    virtual void solveSecondKind(const VectorXd & potential, VectorXd & charge);
    virtual void solveFull(const VectorXd & potential, VectorXd & charge);
    virtual void initPointers();
    element *elementTree; //*E_; Hierarchical element list
    wavelet *waveletList; //*W_; List of wavelets

};
#endif
