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
                GreensFunction * greenInside_;
                GreensFunction * greenOutside_;
		bool allocated;
                virtual std::ostream & printSolver(std::ostream & os) = 0; 
	public:
		PCMSolver() {}
	        PCMSolver(GreensFunction * gfInside, GreensFunction * gfOutside) : greenInside_(gfInside), greenOutside_(gfOutside), allocated(true) {}
                virtual ~PCMSolver()                                                                      
		{
			if (allocated)                                                                    
		        {
				delete greenInside_;                                                       
		       	 	delete greenOutside_;                                                      
		        }                                                                                 
		}                                                                                         
                                                                                                          
                GreensFunction * greenInside() const { return greenInside_; }                              
                GreensFunction * greenOutside() const { return greenOutside_; }                                

		/*! \brief Calculation of the PCM matrix.
		 *  \param[in] cavity the cavity to be used.
		 */
                virtual void buildSystemMatrix(const Cavity & cavity) = 0;
	        /*! \brief Computation of ASC given the MEP. 
		 *  \param[in] potential the vector containing the MEP at cavity points.
		 *  \param[out] charge the vector containing the ASC at cavity points.
		 *  \param[in] irrep the irreducible representation of the MEP and ASC.
		 *
		 *  Given the MEP for a certain irrep, computes the corresponding ASC.
		 *  By default, we expect the totally symmetric irrep to be needed,
		 *  as in energy calculations..
		 */
                virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge, int irrep = 0) = 0; 
                friend std::ostream & operator<<(std::ostream & os, PCMSolver & solver)                      
		{                                                                                         
                    return solver.printSolver(os);                                                           
                }                                                                                         
};

#endif // PCMSOLVER_HPP
