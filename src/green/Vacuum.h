#ifndef VACUUM_H
#define VACUUM_H

#include <iostream>
#include <string>

#include <Eigen/Dense>

#include "GreensFunction.h"
#include "GreensFunctionFactory.h"

/*! \file Vacuum.h
 *  \class Vacuum
 *  \brief Green's functions for vacuum.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013
 */

class Vacuum : public GreensFunction
{
	public:
		Vacuum(int how_) : GreensFunction(how_, true) {} 
		virtual ~Vacuum() {}
 		virtual void compDiagonal(const Eigen::VectorXd & elementArea_, const Eigen::VectorXd & elementRadius_,
                                          Eigen::MatrixXd & S_, Eigen::MatrixXd & D_) const;
	        virtual double getDielectricConstant() const { return 1.0; }	
	private:
		virtual	Eigen::Array4d numericalDirectional(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, 
				 			    Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const;
		virtual	Eigen::Array4d analyticDirectional(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_,
				 			   Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const;
		virtual	Eigen::Array4d automaticDirectional(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_,
							    Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const;  
		virtual	Eigen::Array4d automaticGradient(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, 
				                         Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const;
		virtual	Eigen::Array4d automaticHessian(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, 
				                        Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const;
};

namespace
{
	GreensFunction * createVacuum(int how_, double epsilon_ = 1.0, double kappa_ = 0.0)
	{
		return new Vacuum(how_);
	}
	const std::string VACUUM("Vacuum");
	const bool registeredVacuum = GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(VACUUM, createVacuum);
}

#endif // VACUUM_H
