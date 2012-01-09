/*! \file PCMSolver.h 
\brief PCM solver
*/


#ifndef PCMSOLVER
#define PCMSOLVER

#include <string>
#include <vector>
#include <iostream>
#include <complex>

#include "Solvent.h"

using namespace std;

class PCMSolver{
 public:
    PCMSolver(GreensFunction &gfi, GreensFunction &gfo);
    PCMSolver(GreensFunction *gfi, GreensFunction *gfo);
    PCMSolver(Section solver);
    ~PCMSolver();
    GreensFunction & getGreenInside();
    GreensFunction & getGreenOutside();    
    GreensFunction * getGreenInsideP();
    GreensFunction * getGreenOutsideP();    
    int getCavitySize() const {return cavitySize;};
    void setCavitySize(int size) {cavitySize = size;};

    virtual void buildSystemMatrix(Cavity & cavity) = 0;

    virtual VectorXd compCharge(const VectorXd & potential) = 0;
    virtual void compCharge(const VectorXd & potential, VectorXd & charge) = 0;

    virtual void setSolverType(const string & type);
    virtual void setSolverType(int type);

    virtual int getSolverType() { return solverType; }

    enum solverTypes{Traditional, Wavelet};

    vector<Solvent> initSolvent();
    string & getSolvent(){ return solvent; }
    void setSolvent(string & solvent);

    friend std::ostream & operator<<(std::ostream & o, PCMSolver & c);

 protected:
    bool allocated;
    int cavitySize;
    GreensFunction *greenInside;
    GreensFunction *greenOutside;
    int solverType;
    virtual ostream & printObject(ostream & os);
    string solvent;
};
#endif
