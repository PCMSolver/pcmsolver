#ifndef CPCMSOLVER_H
#define CPCMSOLVER_H

#include <string>
#include <vector>
#include <iostream>
#include <complex>

#include "Config.h"

#include "PCMSolver.h"

/*! \file CPCMSolver.h  
 *  \class CPCMSolver
 *  \brief Solver for conductor-like approximation: C-PCM (COSMO)
 *  \author Roberto Di Remigio
 *  \date 2013
 */

class CPCMSolver : public PCMSolver 
{
	public:
		CPCMSolver(GreensFunctionInterface &gfi, GreensFunctionInterface &gfo);
                CPCMSolver(GreensFunctionInterface *gfi, GreensFunctionInterface *gfo);                
        //        CPCMSolver(const Section & solver);
                ~CPCMSolver();

                const Eigen::MatrixXd &getPCMMatrix() const { return PCMMatrix; }
                                                                                                      
                virtual void buildSystemMatrix(Cavity & cavity);
                virtual void buildAnisotropicMatrix(GePolCavity & cav);
                virtual void buildIsotropicMatrix(GePolCavity & cav);
                //    virtual VectorXd compCharge(const VectorXd & potential);
                virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
                virtual void setCorrection(double _correction);

	private:
    		bool builtAnisotropicMatrix;
    		bool builtIsotropicMatrix;
//    		static const double factor = 1.0694;
    		static const double factor = 1.07;
    		double correction;
    		Eigen::MatrixXd PCMMatrix;
    		virtual std::ostream & printObject(std::ostream & os);
};

#endif // CPCMSOLVER_H
