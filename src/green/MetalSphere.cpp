#include "MetalSphere.hpp"

#include <cmath>
#include <ostream>
#include <stdexcept>

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

extern"C" 
{
void gsfera_cpp_(double* epssol, double* epsre, double* epsim, 
		 double* sphRadius, double* ps, double* p1, double* p2,
		 double* greenre, double* greenim);
}

void MetalSphere::compDiagonal(const Eigen::VectorXd & elementArea_, const Eigen::VectorXd & elementRadius_, Eigen::MatrixXd & S_, Eigen::MatrixXd & D_) const
{
	throw std::runtime_error("Green's function for a metal sphere has not yet been implemented!");
}

Eigen::Array4d MetalSphere::numericalDirectional(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, 
	            			         Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const
{
	Eigen::Array4d result;
	
	// Calculation of the value of the Green's Function
	double epsre, epsim;
	double greenre, greenim;
    	double point1[3], point2[3], sphere[3];
    	point1[0] = source_(0);
    	point1[1] = source_(1);
    	point1[2] = source_(2);
    	point2[0] =  probe_(0);
    	point2[1] =  probe_(1);
    	point2[2] =  probe_(2);
    	sphere[0] = sphPosition(0);
    	sphere[1] = sphPosition(1);
    	sphere[2] = sphPosition(2);
    	epsre = epsMetal.real();
    	epsim = epsMetal.imag();
    	gsfera_cpp_(&epsSolvent, &epsre, &epsim, &sphRadius, sphere, point1, point2, &greenre, &greenim);

	result(0) = greenre;
/* This is for the evaluation of the derivative...
double MetalSphere::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2) {
    return epsSolvent * (this->derivativeProbe(direction, p1, p2));  // NORMALIZTION TEMPORARY REMOVED /direction.norm();
}*/
}

Eigen::Array4d MetalSphere::analyticDirectional(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_, 	
	 			                Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const
{
	throw std::runtime_error("Analytic derivative is not available for MetalSphere.");
}

Eigen::Array4d MetalSphere::automaticDirectional(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_,	
	 			                 Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const
{
	throw std::runtime_error("Automatic directional derivative is not available for MetalSphere.");
}

Eigen::Array4d MetalSphere::automaticGradient(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_,   	
	 			              Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const
{
	throw std::runtime_error("Automatic gradient is not available for MetalSphere.");
}

Eigen::Array4d MetalSphere::automaticHessian(Eigen::Vector3d & sourceNormal_, Eigen::Vector3d & source_,    	
	 			             Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const
{
	throw std::runtime_error("Automatic Hessian is not available for MetalSphere.");
}

std::ostream & MetalSphere::printGreensFunction(std::ostream & os)
{
	return os;
}
