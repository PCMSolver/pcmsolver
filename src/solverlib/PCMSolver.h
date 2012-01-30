/*! \file PCMSolver.h 
\brief PCM solver
*/


#ifndef PCMSOLVER
#define PCMSOLVER

#include <string>
#include <vector>
#include <iostream>
#include <complex>

#include "GreensFunction.h"
#include "Solvent.h"

using namespace std;

template <class T>
class PCMSolver{
 public:
    PCMSolver(GreensFunction<T> &gfi, GreensFunction<T> &gfo);
    PCMSolver(GreensFunction<T> *gfi, GreensFunction<T> *gfo);
    PCMSolver(Section solver);
    ~PCMSolver();
    GreensFunction<T> & getGreenInside();
    GreensFunction<T> & getGreenOutside();    
    GreensFunction<T> * getGreenInsideP();
    GreensFunction<T> * getGreenOutsideP();    
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

    friend std::ostream & operator<<(std::ostream & os, PCMSolver & obj) {
        return obj.printObject(os);
    }

 protected:
    bool allocated;
    int cavitySize;
    GreensFunction<T> *greenInside;
    GreensFunction<T> *greenOutside;
    int solverType;
    virtual ostream & printObject(ostream & os);
    string solvent;
};
#endif
