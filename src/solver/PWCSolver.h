#ifndef PWCSOLVER_H
#define PWCSOLVER_H

#include <iostream>
#include <string>

#include <Eigen/Dense>

#include "Config.h"

class GreensFunction;

extern "C"{
#include "vector3.h"
#include "sparse2.h"
#include "intvector.h"
#include "basis.h"
}

#include "WEMSolver.h"
#include "SolverFactory.h"

/*! \file PWCSolver.h
 *  \class PWCSOlver
 *  \brief Wavelet solver, piecewise constant.
 *  \author Luca Frediani
 *  \date 2012
 */

class PWCSolver : public WEMSolver 
{
	private:
		virtual void initPointers();
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

#endif // PWCSOLVER_H
