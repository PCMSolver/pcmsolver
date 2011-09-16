/*! \file WEMSolver.h 
\brief PCM solver
*/


#ifndef WEMSOLVER
#define WEMSOLVER

#include <string>
#include <vector>
#include <iostream>
#include <complex>

using namespace std;

class WEMSolver : public PCMSolver {
 public:
    WEMSolver(GreensFunction &gfi, GreensFunction &gfo);
    WEMSolver(GreensFunction *gfi, GreensFunction *gfo);
    WEMSolver(Section solver);
    ~WEMSolver();
    virtual VectorXd compCharge(const VectorXd & potential);
    virtual void compCharge(const VectorXd & potential, VectorXd & charge);
};
#endif
