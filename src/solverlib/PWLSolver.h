/*! \file PWLSolver.h 
\brief PCM solver
*/


#ifndef PWLSOLVER_H_
#define PWLSOLVER_H_

#include <string>
#include <vector>
#include <iostream>
#include <complex>

extern "C"{
#include "vector3.h"
#include "sparse2.h"
#include "intvector_pwl.h"
#include "basis_pwl.h"
}

using namespace std;

class PWLSolver : public WEMSolver {
 public:
    PWLSolver(GreensFunctionInterface &gfi, GreensFunctionInterface &gfo);
    PWLSolver(GreensFunctionInterface *gfi, GreensFunctionInterface *gfo);
    PWLSolver(Section solver);
    ~PWLSolver();
 private:
    virtual void initInterpolation();
    virtual void constructWavelets();
    virtual void constructSi();
    virtual void constructSe();
    virtual void solveFirstKind(const VectorXd & potential, VectorXd & charge);
    virtual void solveSecondKind(const VectorXd & potential, VectorXd & charge);
    virtual void solveFull(const VectorXd & potential, VectorXd & charge);
    virtual void initPointers();
    element_pwl *elementTree; //*E_; Hierarchical element list
    wavelet_pwl *waveletList; //*W_; List of wavelets
};
#endif
