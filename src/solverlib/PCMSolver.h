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
    PCMSolver(GreensFunction &gfi, GreensFunction &gfo);
    PCMSolver(GreensFunction *gfi, GreensFunction *gfo);
    PCMSolver(Section solver);
    ~PCMSolver();
    GreensFunction &getGreenInside();
    GreensFunction &getGreenOutside();    
    int getCavitySize() const {return cavitySize;};
    void setCavitySize(int size) {cavitySize = size;};
    virtual VectorXd compCharge(const VectorXd & potential) = 0;
    virtual void compCharge(const VectorXd & potential, VectorXd & charge) = 0;
 protected:
    bool allocated;
    int cavitySize;
    GreensFunction *greenInside;
    GreensFunction *greenOutside;
};
#endif
