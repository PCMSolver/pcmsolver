#ifndef PWCSOLVER_H
#define PWCSOLVER_H

#include <string>
#include <vector>
#include <iostream>
#include <complex>

#include "Config.h"

class GreensFunction;

extern "C"{
#include "vector3.h"
#include "sparse2.h"
#include "intvector.h"
#include "basis.h"
}

#include "WEMSolver.h"

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

#endif
