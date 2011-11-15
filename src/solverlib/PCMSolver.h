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
    const MatrixXd& getPCMMatrix() const {return PCMMatrix;};
    virtual void buildAnisotropicMatrix(GePolCavity cav);
    virtual void buildIsotropicMatrix(GePolCavity cav);
    virtual double compDiagonalElementSoper(GreensFunction *green, int i, 
                                            GePolCavity cav);
    virtual double compDiagonalElementDoper(GreensFunction *green, int i, 
                                            GePolCavity cav);
    virtual VectorXd compCharge(const VectorXd & potential);
    virtual void compCharge(const VectorXd & potential, VectorXd & charge);
 private:
    bool allocated;
    bool builtAnisotropicMatrix;
    bool builtIsotropicMatrix;
    static const double factor = 1.07;
    int cavitySize;
    GreensFunction *greenInside;
    GreensFunction *greenOutside;
    MatrixXd PCMMatrix;
};
#endif
