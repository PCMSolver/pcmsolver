#ifndef CPCMSOLVER_H
#define CPCMSOLVER_H

#include <string>
#include <vector>
#include <iostream>
#include <complex>

#include "Config.h"

class GreensFunction;
class Cavity;
class GePolCavity;

#include "PCMSolver.h"

/*! \file CPCMSolver.h  
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
    		virtual std::ostream & printObject(std::ostream & os);
//    		static const double factor = 1.0694;
    		static const double factor = 1.07;
	public:
		CPCMSolver(GreensFunction &gfi, GreensFunction &gfo, double correction_ = 0.0) 
			: PCMSolver(gfi, gfo), builtIsotropicMatrix(false), builtAnisotropicMatrix(false), correction(correction_) {}
                CPCMSolver(GreensFunction *gfi, GreensFunction *gfo, double correction_ = 0.0) 
			: PCMSolver(gfi, gfo), builtIsotropicMatrix(false), builtAnisotropicMatrix(false), correction(correction_) {}                
                //CPCMSolver(const Section & solver);
                virtual ~CPCMSolver() {}
                const Eigen::MatrixXd & getPCMMatrix() const { return PCMMatrix; }
                virtual void buildSystemMatrix(Cavity & cavity);
                //virtual VectorXd compCharge(const VectorXd & potential);
                virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
                void setCorrection(double correction_) { correction = correction_; }
	private:
                void buildIsotropicMatrix(GePolCavity & cav);
};

#endif // CPCMSOLVER_H
