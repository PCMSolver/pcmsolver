#ifndef CPCMSOLVER_HPP
#define CPCMSOLVER_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

class Cavity;
class GreensFunction;

#include "PCMSolver.hpp"
#include "SolverData.hpp"
#include "SolverFactory.hpp"

/*! \file CPCMSolver.hpp  
 *  \class CPCMSolver
 *  \brief Solver for conductor-like approximation: C-PCM (COSMO)
 *  \author Roberto Di Remigio
 *  \date 2013
 */

class CPCMSolver : public PCMSolver 
{
	private:
    		bool builtIsotropicMatrix;
    		bool builtAnisotropicMatrix;
		bool hermitivitize_;
    		double correction_;
		Eigen::MatrixXd fullPCMMatrix;
		std::vector<Eigen::MatrixXd> blockPCMMatrix;
                void buildIsotropicMatrix(const Cavity & cavity);
    		virtual std::ostream & printSolver(std::ostream & os);
	public:
		CPCMSolver() {}
                CPCMSolver(GreensFunction * gfInside, GreensFunction * gfOutside, bool symm, double correction) 
			: PCMSolver(gfInside, gfOutside), builtIsotropicMatrix(false), builtAnisotropicMatrix(false), hermitivitize_(symm), correction_(correction) {}                
                virtual ~CPCMSolver() {}
                virtual void buildSystemMatrix(const Cavity & cavity);
                virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge, int irrep = 0);
                void correction(double corr) { correction_ = corr; }
                friend std::ostream & operator<<(std::ostream & os, CPCMSolver & solver)                      
		{                                                                                         
                    return solver.printSolver(os);                                                           
                }                                                                                         
};

namespace
{
	PCMSolver * createCPCMSolver(const solverData & _data)
	{
		return new CPCMSolver(_data.gfInside, _data.gfOutside, _data.hermitivitize, _data.correction);
	}
	const std::string CPCMSOLVER("CPCM");
	const bool registeredCPCMSolver = SolverFactory::TheSolverFactory().registerSolver(CPCMSOLVER, createCPCMSolver);
}

#endif // CPCMSOLVER_HPP
