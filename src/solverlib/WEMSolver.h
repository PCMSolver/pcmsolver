/*! \file WEMSolver.h 
\brief PCM solver
*/


#ifndef WEMSOLVER
#define WEMSOLVER

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

class WEMSolver : public PCMSolver {
 public:
    WEMSolver(GreensFunction &gfi, GreensFunction &gfo);
    WEMSolver(GreensFunction *gfi, GreensFunction *gfo);
    WEMSolver(Section solver);
    ~WEMSolver();
    virtual void constructSystemMatrix();
    virtual VectorXd compCharge(const VectorXd & potential);
    virtual void compCharge(const VectorXd & potential, VectorXd & charge);
 protected:
  // Parameters
  double eps_;
  unsigned int quadratureLevel_;
  // System matrices
  sparse2 S_i_, S_e_;
  bool systemMatricesInitialized_;
  // Point list
  vector3 *P_;
  // Element list
  unsigned int **F_;
  // Something 1
  vector3 ****T_;
  // Number of knot points or something
  unsigned int np_;
  // Number of ansatz functions
  unsigned int nf_;
  // Number of points 
  unsigned int p_;
  // Patch level (2**M * 2**M elements per patch)
  unsigned int M_;
  // Hierarchical element list
  element *E_;
  // List of wavelets
  wavelet *W_;
  // Number of quadrature points
  int nPoints_;

  // System matrix sparsities
  double apriori1_, aposteriori1_;
  double apriori2_, aposteriori2_;
};
#endif
