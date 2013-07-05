#ifndef PCMSOLVER_HPP
#define PCMSOLVER_HPP

#include <iosfwd>

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
	protected:
                GreensFunction * greenInside;
                GreensFunction * greenOutside;
		bool allocated;
                virtual std::ostream & printSolver(std::ostream & os) = 0; 
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
                                                                                                          
                GreensFunction * getGreenInside() const { return greenInside; }                              
                GreensFunction * getGreenOutside() const { return greenOutside; }                                
                                                                                                          
                virtual void buildSystemMatrix(Cavity & cavity) = 0;                                      
                // ask jonas virtual VectorXd compCharge(const VectorXd & potential);                     
                virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge) = 0; 
                friend std::ostream & operator<<(std::ostream & os, PCMSolver & solver)                      
		{                                                                                         
                    return solver.printSolver(os);                                                           
                }                                                                                         
};

#endif // PCMSOLVER_HPP
