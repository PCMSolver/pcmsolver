#ifndef PCMSOLVER_H
#define PCMSOLVER_H

#include <string>
#include <vector>
#include <iostream>
#include <complex>

#include "Config.h"

#include "GreensFunction.h"
#include "Solvent.h"

/*! 
 * \file PCMSolver.h
 * \class PCMSolver
 * \brief Abstract Base Class for solvers inheritance hierarchy.
 * \author Luca Frediani 
 * \date 2011
 */

class PCMSolver
{
	public:
		 PCMSolver(GreensFunctionInterface &gfi, GreensFunctionInterface &gfo, 
              		int equation=SecondKind, int solver=IEFPCM);
	         PCMSolver(GreensFunctionInterface *gfi, GreensFunctionInterface *gfo, 
             	 	int equation=SecondKind, int solver=IEFPCM);
//                 PCMSolver(const Section & solver);                                                        
                 ~PCMSolver();

                 GreensFunctionInterface & getGreenInside();
                 GreensFunctionInterface & getGreenOutside();    
                 GreensFunctionInterface * getGreenInsideP();
                 GreensFunctionInterface * getGreenOutsideP();    

                 int getCavitySize() const {return cavitySize;};
                 void setCavitySize(int size) {cavitySize = size;};
                 virtual void buildSystemMatrix(Cavity & cavity) = 0;
                 // ask jonas virtual VectorXd compCharge(const VectorXd & potential);
                 virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge) = 0;
                 virtual void setSolverType(const std::string & type);
                 virtual void setSolverType(int type);
                 virtual int getSolverType() { return solverType; }
                 virtual void setEquationType(const std::string & type);
                 virtual void setEquationType(int type);
                 virtual int getEquationType() { return integralEquation; }
                 Solvent getSolvent() { return solvent; }
                 void setSolvent(Solvent & solvent);
                 friend std::ostream & operator<<(std::ostream & os, PCMSolver & obj) 
		 {
                     return obj.printObject(os);
                 }
                 enum EquationType {FirstKind, SecondKind, Full};
                 enum SolverType {IEFPCM, CPCM, Wavelet, Linear};
                 bool isPWL(){return (solverType == Linear);}
	
	protected:
		 bool allocated;
                 int cavitySize;                              
                 GreensFunctionInterface *greenInside;
                 GreensFunctionInterface *greenOutside;
                 SolverType solverType;
                 EquationType integralEquation;
                 virtual std::ostream & printObject(std::ostream & os);
                 Solvent solvent;
};

#endif // PCMSOLVER_H
