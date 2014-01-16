#ifndef CPCMSOLVER_HPP
#define CPCMSOLVER_HPP

#include <iosfwd>
#include <string>

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
    		double correction;
    		Eigen::MatrixXd PCMMatrix;
//    		static const double factor = 1.0694;
//    		static const double factor = 1.07;
                void buildIsotropicMatrix(Cavity & cav);
    		virtual std::ostream & printSolver(std::ostream & os);
	public:
		CPCMSolver() {}
                CPCMSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_, double correction_ = 0.0) 
			: PCMSolver(gfInside_, gfOutside_), builtIsotropicMatrix(false), builtAnisotropicMatrix(false), correction(correction_) {}                
                //CPCMSolver(const Section & solver);
                virtual ~CPCMSolver() {}
                const Eigen::MatrixXd & getPCMMatrix() const { return PCMMatrix; }
                virtual void buildSystemMatrix(Cavity & cavity);
                //virtual VectorXd compCharge(const VectorXd & potential);
                virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
                void setCorrection(double correction_) { correction = correction_; }
                friend std::ostream & operator<<(std::ostream & os, CPCMSolver & solver)                      
		{                                                                                         
                    return solver.printSolver(os);                                                           
                }                                                                                         
};

namespace
{
	PCMSolver * createCPCMSolver(const solverData & _data)
	{
		return new CPCMSolver(_data.gfInside, _data.gfOutside, _data.correction);
	}
	const std::string CPCMSOLVER("CPCM");
	const bool registeredCPCMSolver = SolverFactory::TheSolverFactory().registerSolver(CPCMSOLVER, createCPCMSolver);
}

#endif // CPCMSOLVER_HPP
