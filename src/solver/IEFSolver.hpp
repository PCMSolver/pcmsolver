#ifndef IEFSOLVER_HPP
#define IEFSOLVER_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

class Cavity;
class GreensFunction;

#include "PCMSolver.hpp"
#include "SolverData.hpp"
#include "SolverFactory.hpp"

/*! 
 * \file IEFSolver.hpp
 * \class IEFSolver
 * \brief Traditional solver.
 * \author Luca Frediani 
 * \date 2011
 */

class IEFSolver : public PCMSolver 
{
	private:
   	 	bool builtIsotropicMatrix;
    		bool builtAnisotropicMatrix;
//	 	static const double factor = 1.0694;
//    		static const double factor = 1.07;
		Eigen::MatrixXd PCMMatrix;
                void buildAnisotropicMatrix(Cavity & cav);
                void buildIsotropicMatrix(Cavity & cav);
    		virtual std::ostream & printSolver(std::ostream & os);
	public:
		IEFSolver() {}
    		IEFSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_) 
			: PCMSolver(gfInside_, gfOutside_), builtIsotropicMatrix(false), builtAnisotropicMatrix(false) {}
                virtual ~IEFSolver() {}
                const Eigen::MatrixXd & getPCMMatrix() const { return PCMMatrix; }
                virtual void buildSystemMatrix(Cavity & cavity);
                //virtual VectorXd compCharge(const VectorXd & potential);
                virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
                friend std::ostream & operator<<(std::ostream & os, IEFSolver & solver)                      
		{                                                                                         
                    return solver.printSolver(os);                                                           
                }                                                                                         
};

namespace
{
	PCMSolver * createIEFSolver(const solverData & _data)
	{
		return new IEFSolver(_data.gfInside, _data.gfOutside);
	}
	const std::string IEFSOLVER("IEFPCM");
	const bool registeredIEFSolver = SolverFactory::TheSolverFactory().registerSolver(IEFSOLVER, createIEFSolver);
}

#endif // IEFSOLVER_HPP
