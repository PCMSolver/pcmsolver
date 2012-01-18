/*! \file WEMSolver.h 
\brief PCM solver
*/


#ifndef WEMSOLVER_H_
#define WEMSOLVER_H_

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
    vector3 **** getT_(){return T_;}
    int getQuadratureLevel(){return quadratureLevel_;}

    virtual void buildSystemMatrix(Cavity & cavity);
    
    virtual void constructSystemMatrix();
    virtual void uploadCavity(WaveletCavity cavity);
    virtual VectorXd compCharge(const VectorXd & potential);
    virtual void compCharge(const VectorXd & potential, VectorXd & charge);
 protected:
  // Parameters
    
    double threshold;
    unsigned int quadratureLevel_;
    // System matrices
    sparse2 S_i_, S_e_;
    bool systemMatricesInitialized_;
    // Point list
    vector3 *nodeList; //*P_; --
    // Element list
    unsigned int **elementList; //**F_; --
    // Something 1 
    vector3 ****T_; //--
    // Number of knot points or something
    unsigned int nNodes; //np_; --
    // Number of ansatz functions
    unsigned int nFunctions; //nf_; --
    // Number of points 
    unsigned int nPatches; // p_; --
    // Patch level (2**M * 2**M elements per patch)
    unsigned int nLevels; //M_; --
    // Hierarchical element list
    element *elementTree; //*E_;
    // List of wavelets
    wavelet *waveletList; //*W_;
    // Number of quadrature points
    int nQuadPoints; // nPoints_;
    
    // System matrix sparsities
    double apriori1_, aposteriori1_;
    double apriori2_, aposteriori2_;
};
#endif
