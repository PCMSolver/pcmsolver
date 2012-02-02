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
    PWCSolver(Section solver);
    ~PWCSolver();
    virtual void constructSystemMatrix();
    virtual void initInterpolation();
    virtual void initPointers();
    virtual void compCharge(const VectorXd & potential, VectorXd & charge);
 protected:
    element *elementTree; //*E_; Hierarchical element list
    wavelet *waveletList; //*W_; List of wavelets
};
#endif
