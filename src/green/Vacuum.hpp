#ifndef VACUUM_HPP
#define VACUUM_HPP

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

#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file Vacuum.hpp
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
	        virtual double dielectricConstant() const { return 1.0; }	
                friend std::ostream & operator<<(std::ostream & os, Vacuum & green)                      
		{                                                                                         
                    return green.printGreensFunction(os);                                                           
                }                                                                                         
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
                virtual std::ostream & printGreensFunction(std::ostream & os); 
};

namespace
{
	GreensFunction * createVacuum(const greenData & _data)
	{
		return new Vacuum(_data.how);
	}
	const std::string VACUUM("Vacuum");
	const bool registeredVacuum = GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(VACUUM, createVacuum);
}

#endif // VACUUM_HPP
