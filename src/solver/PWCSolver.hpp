#ifndef PWCSOLVER_HPP
#define PWCSOLVER_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

extern "C"
{
#include "vector3.h"
#include "sparse2.h"
#include "intvector.h"
#include "basis.h"
}

class GreensFunction;

#include "SolverFactory.hpp"
#include "WEMSolver.hpp"

/*! \file PWCSolver.hpp
 *  \class PWCSOlver
 *  \brief Wavelet solver, piecewise constant.
 *  \author Luca Frediani
 *  \date 2012
 */

class PWCSolver : public WEMSolver 
{
	private:
		virtual void initPointers();
    		virtual std::ostream & printSolver(std::ostream & os);
	public:
                PWCSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_) : WEMSolver(gfInside_, gfOutside_, Full)
		{
	                initPointers();
		}
                PWCSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_, int integralEquation_) : WEMSolver(gfInside_, gfOutside_, integralEquation_)
		{
	                initPointers();
		}
		//PWCSolver(const Section & solver);
                virtual ~PWCSolver();
                friend std::ostream & operator<<(std::ostream & os, PWCSolver & solver)                      
		{                                                                                         
                    return solver.printSolver(os);                                                           
                }                                                                                         
	private: 
		virtual void initInterpolation();
                virtual void constructWavelets();                                            
                virtual void constructSi();
                virtual void constructSe();
                virtual void solveFirstKind(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
                virtual void solveSecondKind(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
                virtual void solveFull(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
                element *elementTree; //*E_; Hierarchical element list
                wavelet *waveletList; //*W_; List of wavelets
};

namespace
{
	PCMSolver * createPWCSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_, double correction_ = 0.0, int integralEquation_ = 1)
	{
		return new PWCSolver(gfInside_, gfOutside_, integralEquation_);
	}
	const std::string PWCSOLVER("Wavelet");
	const bool registeredPWCSolver = SolverFactory::TheSolverFactory().registerSolver(PWCSOLVER, createPWCSolver);
}

#endif // PWCSOLVER_HPP
