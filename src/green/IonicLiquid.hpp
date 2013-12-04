#ifndef IONICLIQUID_HPP
#define IONICLIQUID_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file IonicLiquid.hpp
 *  \class IonicLiquid
 *  \brief Green's functions for ionic liquid, described by the linearized Poisson-Boltzmann equation.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013
 */

class IonicLiquid : public GreensFunction
{
	public:
		IonicLiquid(int how_, double epsilon_, double kappa_) : GreensFunction(how_), epsilon(epsilon_), kappa(kappa_) {}
		virtual ~IonicLiquid() {}
 		virtual void compDiagonal(const Eigen::VectorXd & elementArea_, const Eigen::VectorXd & elementRadius_,
                                          Eigen::MatrixXd & S_, Eigen::MatrixXd & D_) const; 
	       	virtual double getDielectricConstant() const { return epsilon; }
                friend std::ostream & operator<<(std::ostream & os, IonicLiquid & green)                      
		{                                                                                         
                    return green.printGreensFunction(os);                                                           
                }                                                                                         
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
                virtual std::ostream & printGreensFunction(std::ostream & os); 
};

namespace
{
	GreensFunction * createIonicLiquid(int how_, double epsilon_ = 1.0, double kappa_ = 0.0)
	{
		return new IonicLiquid(how_, epsilon_, kappa_);
	}
	const std::string IONICLIQUID("IonicLiquid");
	const bool registeredIonicLiquid = GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(IONICLIQUID, createIonicLiquid);
}

#endif // IONICLIQUID_HPP
