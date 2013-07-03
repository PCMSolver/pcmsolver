#ifndef PCMSOLVER_H
#define PCMSOLVER_H

#include <string>
#include <vector>
#include <iostream>
#include <complex>

#include "Config.h"

class Cavity;

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
		PCMSolver() {}
		PCMSolver(GreensFunction &gfi, GreensFunction &gfo, int equation = SecondKind, int solver = IEFPCM)
			: greenInside(&gfi), greenOutside(&gfo), solverType(solver), integralEquation(equation), allocated(false) {}
	        PCMSolver(GreensFunction *gfi, GreensFunction *gfo, int equation = SecondKind, int solver = IEFPCM)
			: greenInside(gfi), greenOutside(gfo), solverType(solver), integralEquation(equation), allocated(false) {}
//                 PCMSolver(const Section & solver);                                                       
                virtual ~PCMSolver()                                                                      
		{                                                                                         
		        if (allocated)                                                                    
		        {                                                                                 
		       	 delete greenInside;                                                       
		       	 delete greenOutside;                                                      
		        }                                                                                 
		}                                                                                         
                                                                                                          
                GreensFunction & getGreenInside() { return *greenInside; }			
                GreensFunction & getGreenOutside() { return *greenOutside; }
                GreensFunction * getGreenInsideP() { return greenInside; }                              
                GreensFunction * getGreenOutsideP() { return greenOutside; }                                
                                                                                                          
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
                void setSolvent(Solvent & solvent_) { solvent = solvent_; }
                friend std::ostream & operator<<(std::ostream & os, PCMSolver & obj)                      
		{                                                                                         
                    return obj.printObject(os);                                                           
                }                                                                                         
                bool isPWL(){return (solverType == Linear);}                                              
	
	protected:
                GreensFunction *greenInside;
                GreensFunction *greenOutside;
                int solverType;
                int integralEquation;
		bool allocated;
                virtual std::ostream & printObject(std::ostream & os);
                Solvent solvent;
                enum EquationType {FirstKind, SecondKind, Full};                                          
                enum SolverType {IEFPCM, CPCM, Wavelet, Linear};                                          
};

#endif // PCMSOLVER_H
