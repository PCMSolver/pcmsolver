#ifndef METALSPHERE_HPP
#define METALSPHERE_HPP

#include <iosfwd>
#include <complex>

#include "Config.hpp"

#include <Eigen/Dense>

#include "GreensFunction.hpp"

/*! \file MetalSphere.hpp
 *  \class MetalSphere
 *  \brief Class to describe spherical metal nanoparticles.
 *  \author Luca Frediani
 *  \date 2011
 *
 *  This class is a wrapper around the FORTRAN routines written
 *  by Stefano Corni et al. to take into account the presence
 *  of metal nanoparticles.
 */


class MetalSphere : public GreensFunction
{
	private:
		typedef std::complex<double> dcomplex;
	public:
		MetalSphere(int how_, double eps_, double epsRe_, double epsIm_, Eigen::Vector3d & pos_, double radius_)
			: GreensFunction(how_, false), epsSolvent(eps_), epsMetal(dcomplex(epsRe_, epsIm_)), sphPosition(pos_), sphRadius(radius_) 
                virtual ~MetalSphere() {}                                              
 		virtual void compDiagonal(const Eigen::VectorXd & elementArea_, const Eigen::VectorXd & elementRadius_,
                                          Eigen::MatrixXd & S_, Eigen::MatrixXd & D_) const;
	        virtual double getDielectricConstant() const {}	
                friend std::ostream & operator<<(std::ostream & os, MetalSphere & green)                      
		{                                                                                         
                    return green.printGreensFunction(os);                                                           
                }                                                                                         
	private:
                double epsSolvent;
                dcomplex epsMetal;
                Eigen::Vector3d sphPosition;
                double sphRadius;
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
#endif // METALSPHERE_HPP
