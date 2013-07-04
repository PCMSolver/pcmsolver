#ifndef UNIFORMDIELECTRIC_HPP
#define UNIFORMDIELECTRIC_HPP

#include <iostream>
#include <string>

#include <Eigen/Dense>

#include "Config.hpp"

#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file UniformDielectric.hpp
 *  \class UniformDielectric
 *  \brief Green's functions for uniform dielectric.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013
 */

class UniformDielectric : public GreensFunction
{
	public:
		UniformDielectric(int how_, double epsilon_) : GreensFunction(how_, true), epsilon(epsilon_) {}
		virtual ~UniformDielectric() {}
 		virtual void compDiagonal(const Eigen::VectorXd & elementArea_, const Eigen::VectorXd & elementRadius_,
                                          Eigen::MatrixXd & S_, Eigen::MatrixXd & D_) const;
	       	virtual double getDielectricConstant() const { return epsilon; }
	private:
		double epsilon;
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
	GreensFunction * createUniformDielectric(int how_, double epsilon_ = 1.0, double kappa_ = 0.0)
	{
		return new UniformDielectric(how_, epsilon_);
	}
	const std::string UNIFORMDIELECTRIC("UniformDielectric");
	const bool registeredUniformDielectric = GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(UNIFORMDIELECTRIC, createUniformDielectric);
}

#endif // UNIFORMDIELECTRIC_HPP
