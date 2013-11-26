#ifndef GREENSFUNCTION_HPP
#define GREENSFUNCTION_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

/*! \file GreensFunction.hpp
 *  \class GreensFunction
 *  \brief An Abstract Base Class for Green's functions
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013
 */

class GreensFunction
{
	public:
		GreensFunction(int how_) : how(how_), uniform(false) {}
		GreensFunction(int how_, bool uniform_) : how(how_), uniform(uniform_) {}
		virtual ~GreensFunction() {}
		Eigen::Array4d evaluate(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const; 
		/*! \brief Compute the off-diagonal elements of the S and D matrices.
                 *  \param[in] elementCenter_ the matrix containing the centers of the finite elements.
                 *  \param[in] elementNormal_ the matrix containing the normal vectors to the centers of the finite elements.
                 *  \param[out] S_ the S matrix.
                 *  \param[out] D_ the D matrix.
		 */
		void compOffDiagonal(const Eigen::Matrix3Xd & elementCenter_, const Eigen::Matrix3Xd & elementNormal_, 
                                             Eigen::MatrixXd & S_, Eigen::MatrixXd & D_) const;
		/*! \brief Compute diagonal elements of the S and D matrices.
                 *  \param[in] elementArea_ the vector containing the areas of the finite elements.
                 *  \param[in] elementRadius_ the vector containing the radii of the spheres the finite element belongs to.
                 *  \param[out] S_ the S matrix.
                 *  \param[out] D_ the D matrix.
                 */
 		virtual void compDiagonal(const Eigen::VectorXd & elementArea_, const Eigen::VectorXd & elementRadius_,
                                          Eigen::MatrixXd & S_, Eigen::MatrixXd & D_) const = 0;
		bool isUniform() const { return uniform; }
		virtual double getDielectricConstant() const = 0;
                friend std::ostream & operator<<(std::ostream & os, GreensFunction & green)                      
		{                                                                                         
                    return green.printGreensFunction(os);                                                           
                }                                                                                         
	protected:
		int how;
		bool uniform;
  	private:		
	        /*! \brief Numerical evaluation strategy.
		 *  \param[in] sourceNormal_ the normal vector relative to the source point.
		 *  \param[in] source_ the source point.
		 *  \param[in] probeNormal_ the normal vector relative to the probe point.
		 *  \param[in] probe_ the probe point.
		 */
		virtual Eigen::Array4d numericalDirectional(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, 
							    Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const = 0;
	        /*! \brief Analytic evaluation strategy.
		 *  \param[in] sourceNormal_ the normal vector relative to the source point.
		 *  \param[in] source_ the source point.
		 *  \param[in] probeNormal_ the normal vector relative to the probe point.
		 *  \param[in] probe_ the probe point.
		 */
		virtual Eigen::Array4d analyticDirectional(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, 
         						   Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const = 0;
	        /*! \brief Automatic differentiation evaluation strategy. Calculates directional derivatives internally.
		 *  \param[in] sourceNormal_ the normal vector relative to the source point.
		 *  \param[in] source_ the source point.
		 *  \param[in] probeNormal_ the normal vector relative to the probe point.
		 *  \param[in] probe_ the probe point.
		 */
		virtual Eigen::Array4d automaticDirectional(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, 
							    Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const = 0;
	        /*! \brief Automatic differentiation evaluation strategy. Calculates the full gradient internally.
		 *  \param[in] sourceNormal_ the normal vector relative to the source point.
		 *  \param[in] source_ the source point.
		 *  \param[in] probeNormal_ the normal vector relative to the probe point.
		 *  \param[in] probe_ the probe point.
		 */
		virtual Eigen::Array4d automaticGradient(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, 
						         Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const = 0;
	        /*! \brief Automatic differentiation evaluation strategy. Calculates the full gradient and Hessian internally.
		 *  \param[in] sourceNormal_ the normal vector relative to the source point.
		 *  \param[in] source_ the source point.
		 *  \param[in] probeNormal_ the normal vector relative to the probe point.
		 *  \param[in] probe_ the probe point.
		 */
		virtual Eigen::Array4d automaticHessian(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, 
							Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const = 0;
                virtual std::ostream & printGreensFunction(std::ostream & os) = 0; 
};

#endif // GREENSFUNCTION_HPP
