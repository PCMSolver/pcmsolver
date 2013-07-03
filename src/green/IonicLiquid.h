#ifndef IONICLIQUID_H
#define IONICLIQUID_H

#include <iostream>
#include <string>

#include <Eigen/Dense>

#include "GreensFunction.h"
#include "GreensFunctionFactory.h"

/*! \file IonicLiquid.h
 *  \class IonicLiquid
 *  \brief Green's functions for ionic liquid, described by the linearized Poisson-Boltzmann equation.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013
 */

class IonicLiquid : public GreensFunction
{
	public:
		IonicLiquid(const std::string & how_, double epsilon_, double kappa_) : GreensFunction(how_), epsilon(epsilon_), kappa(kappa_) {}
		virtual ~IonicLiquid() {}
 		virtual void compDiagonal(const Eigen::VectorXd & elementArea_, const Eigen::VectorXd & elementRadius_,
                                          Eigen::VectorXd & diagonalS_, Eigen::VectorXd & diagonalD_) const; 
	private:
		double epsilon;
		double kappa;
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
	GreensFunction * createIonicLiquid(const std::string & how_, double epsilon_ = 1.0, double kappa_ = 0.0)
	{
		return new IonicLiquid(how_, epsilon_, kappa_);
	}
	const std::string IONICLIQUID("IonicLiquid");
	const bool registeredIonicLiquid = GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(IONICLIQUID, createIonicLiquid);
}

#endif // IONICLIQUID_H
