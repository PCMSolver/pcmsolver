#ifndef IEFSOLVER_H
#define IEFSOLVER_H

#include <string>
#include <vector>
#include <iostream>
#include <complex>

/*! 
 * \file IEFSolver.h
 * \class IEFSolver
 * \brief Traditional solver.
 * \author Luca Frediani 
 */

class IEFSolver : public PCMSolver 
{
	public:
		IEFSolver(GreensFunctionInterface &gfi, GreensFunctionInterface &gfo);
    		IEFSolver(GreensFunctionInterface *gfi, GreensFunctionInterface *gfo);
                IEFSolver(const Section & solver);                                      
                ~IEFSolver();

                const Eigen::MatrixXd & getPCMMatrix() const {return PCMMatrix;};
                                                                                        
                virtual void buildSystemMatrix(Cavity & cavity);
                virtual void buildAnisotropicMatrix(GePolCavity & cav);
                virtual void buildIsotropicMatrix(GePolCavity & cav);
                //    virtual VectorXd compCharge(const VectorXd & potential);
                virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
	
	private:
    		bool builtAnisotropicMatrix;
   	 	bool builtIsotropicMatrix;
//	 	static const double factor = 1.0694;
    		static const double factor = 1.07;
		Eigen::MatrixXd PCMMatrix;
    		virtual std::ostream & printObject(std::ostream & os);
};
#endif
