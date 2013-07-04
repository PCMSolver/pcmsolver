#ifndef PWLSOLVER_H
#define PWLSOLVER_H

#include <string>
#include <vector>
#include <iostream>
#include <complex>

#include "Config.h"

class GreensFunction;

extern "C"{
#include "vector3.h"
#include "sparse2.h"
#include "intvector_pwl.h"
#include "basis_pwl.h"
}

#include "WEMSolver.h"

/*! \file PWLSolver.h
 *  \class PWLSOlver
 *  \brief Wavelet solver, piecewise linear.
 *  \author Luca Frediani
 *  \date 2012
 */

class PWLSolver : public WEMSolver 
{
	private:
		virtual void initPointers();
	public:
                PWLSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_) : WEMSolver(gfInside_, gfOutside_, FirstKind)
		{
	                initPointers();
		}
                PWLSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_, int integralEquation_) : WEMSolver(gfInside_, gfOutside_, integralEquation_)
		{
	                initPointers();
		}
		//PWLSolver(Section solver);
                virtual ~PWLSolver();

 	private:
		virtual void initInterpolation();
                virtual void constructWavelets();                                            
                virtual void constructSi();
                virtual void constructSe();
                virtual void solveFirstKind(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
                virtual void solveSecondKind(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
                virtual void solveFull(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
                element_pwl *elementTree; //*E_; Hierarchical element list
                wavelet_pwl *waveletList; //*W_; List of wavelets
};
#endif
