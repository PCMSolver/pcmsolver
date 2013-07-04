#ifndef IEFSOLVER_H
#define IEFSOLVER_H

#include <string>
#include <vector>
#include <iostream>
#include <complex>

#include "Config.h"

class GreensFunction;
class Cavity;
class GePolCavity;

#include "PCMSolver.h"
#include "SolverFactory.h"

/*! 
 * \file IEFSolver.h
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
    		static const double factor = 1.07;
		Eigen::MatrixXd PCMMatrix;
    		virtual std::ostream & printObject(std::ostream & os);
	public:
    		IEFSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_) 
			: PCMSolver(gfInside_, gfOutside_), builtIsotropicMatrix(false), builtAnisotropicMatrix(false) {}
                virtual ~IEFSolver() {}
                const Eigen::MatrixXd & getPCMMatrix() const { return PCMMatrix; }
                virtual void buildSystemMatrix(Cavity & cavity);
                //virtual VectorXd compCharge(const VectorXd & potential);
                virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
	private:
                void buildAnisotropicMatrix(GePolCavity & cav);
                void buildIsotropicMatrix(GePolCavity & cav);
	
};

namespace
{
	PCMSolver * createIEFSolver(GreensFunction * gfInside_, GreensFunction * gfOutside_, double correction_ = 0.0, int integralEquation_ = 1)
	{
		return new IEFSolver(gfInside_, gfOutside_);
	}
	const std::string IEFSOLVER("IEFPCM");
	const bool registeredIEFSolver = SolverFactory::TheSolverFactory().registerSolver(IEFSOLVER, createIEFSolver);
}

#endif // IEFSOLVER_H
