#ifndef PCMSOLVER_H
#define PCMSOLVER_H

#include <string>
#include <vector>
#include <iostream>
#include <complex>

#include "Config.h"

class Cavity;

#include "GreensFunction.h"

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
		PCMSolver() {}
	        PCMSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_)//, int solverType_ = IEFPCM)
			: greenInside(gfInside_), greenOutside(gfOutside_), allocated(true) {} //, solverType(solverType_), allocated(false) {}
                virtual ~PCMSolver()                                                                      
		{
			if (allocated)                                                                    
		        {
				delete greenInside;                                                       
		       	 	delete greenOutside;                                                      
		        }                                                                                 
		}                                                                                         
                                                                                                          
                GreensFunction * getGreenInside() { return greenInside; }                              
                GreensFunction * getGreenOutside() { return greenOutside; }                                
                                                                                                          
                virtual void buildSystemMatrix(Cavity & cavity) = 0;                                      
                // ask jonas virtual VectorXd compCharge(const VectorXd & potential);                     
                virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge) = 0; 
               /* virtual void setSolverType(const std::string & type);                                     
                virtual void setSolverType(int type);                                                     
                virtual int getSolverType() { return solverType; }                                        
                virtual void setEquationType(const std::string & type);                                   
                virtual void setEquationType(int type);                                                   
                virtual int getEquationType() { return integralEquation; } */
                friend std::ostream & operator<<(std::ostream & os, PCMSolver & obj)                      
		{                                                                                         
                    return obj.printObject(os);                                                           
                }                                                                                         
                //bool isPWL(){return (solverType == Linear);}                                              
	protected:
                GreensFunction * greenInside;
                GreensFunction * greenOutside;
                //int solverType;
		bool allocated;
                virtual std::ostream & printObject(std::ostream & os) {}
                //enum SolverType {IEFPCM, CPCM, Wavelet, Linear};                                          
};

#endif // PCMSOLVER_H
