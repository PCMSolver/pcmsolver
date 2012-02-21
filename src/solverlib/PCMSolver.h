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

class PCMSolver{
 public:
    PCMSolver(GreensFunctionInterface &gfi, GreensFunctionInterface &gfo, 
              int equation=SecondKind, int solver=Traditional);
    PCMSolver(GreensFunctionInterface *gfi, GreensFunctionInterface *gfo, 
              int equation=SecondKind, int solver=Traditional);
    PCMSolver(const Section & solver);
    ~PCMSolver();
    GreensFunctionInterface & getGreenInside();
    GreensFunctionInterface & getGreenOutside();    
    GreensFunctionInterface * getGreenInsideP();
    GreensFunctionInterface * getGreenOutsideP();    
    int getCavitySize() const {return cavitySize;};
    void setCavitySize(int size) {cavitySize = size;};
    virtual void buildSystemMatrix(Cavity & cavity) = 0;
    // ask jonas virtual VectorXd compCharge(const VectorXd & potential);
    virtual void compCharge(const VectorXd & potential, VectorXd & charge) = 0;
    virtual void setSolverType(const string & type);
    virtual void setSolverType(int type);
    virtual int getSolverType() { return solverType; }
    virtual void setEquationType(const string & type);
    virtual void setEquationType(int type);
    virtual int getEquationType() { return integralEquation; }
    vector<Solvent> initSolvent();
    Solvent getSolvent(){ return solvent; }
    void setSolvent(Solvent & solvent);
    friend std::ostream & operator<<(std::ostream & os, PCMSolver & obj) {
        return obj.printObject(os);
    }
    enum EquationType {FirstKind, SecondKind, Full};
    enum SolverType {Traditional, Wavelet, Linear};
    bool isPWL(){return (solverType == Linear);}
 protected:
    bool allocated;
    int cavitySize;
    GreensFunctionInterface *greenInside;
    GreensFunctionInterface *greenOutside;
    SolverType solverType;
    EquationType integralEquation;
    virtual ostream & printObject(ostream & os);
    Solvent solvent;
};
#endif
