#ifndef PCMSOLVER_HPP
#define PCMSOLVER_HPP

#include <iostream>

#include "Config.hpp"

class Cavity;

#include "GreensFunction.hpp"

/*! 
 * \file PCMSolver.hpp
 * \class PCMSolver
 * \brief Abstract Base Class for solvers inheritance hierarchy.
 * \author Luca Frediani 
 * \date 2011
 */

class PCMSolver
{
	public:
		PCMSolver() {}
	        PCMSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_) : greenInside(gfInside_), greenOutside(gfOutside_), allocated(true) {}
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
                friend std::ostream & operator<<(std::ostream & os, PCMSolver & obj)                      
		{                                                                                         
                    return obj.printObject(os);                                                           
                }                                                                                         
	protected:
                GreensFunction * greenInside;
                GreensFunction * greenOutside;
		bool allocated;
                virtual std::ostream & printObject(std::ostream & os) {}
};

#endif // PCMSOLVER_HPP
