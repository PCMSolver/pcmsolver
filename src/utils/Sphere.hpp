#ifndef SPHERE_HPP
#define SPHERE_HPP

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

/*! \file Sphere.hpp
 *  \class Sphere
 *  \brief Class describing a sphere.
 *  \author Roberto Di Remigio
 *  \date 2011
 */

class Sphere
{
  	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW // See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
		Sphere() {}
		Sphere(const Eigen::Vector3d & center, double radius, const std::string & colour = "Violet" )
			: sphereCenter_(center), sphereRadius_(radius), sphereColour_(colour) {}
		~Sphere() {}
                double sphereRadius() const { return sphereRadius_; }
                void sphereRadius(double radius) { sphereRadius_ = radius; }
		const Eigen::Vector3d & sphereCenter() const { return sphereCenter_; }
                double sphereCenter(int i) const { return sphereCenter_(i); }
                void sphereCenter(Eigen::Vector3d & coord) { sphereCenter_ = coord; }
		const std::string & sphereColour() const { return sphereColour_; }
                void sphereColour(std::string & colour) { sphereColour_ = colour; }
         	friend inline void swap(Sphere & left, Sphere & right);
         	inline void swap(Sphere & other);
	        /// Assignment operator.
                Sphere & operator=(Sphere other);
                friend std::ostream& operator<<(std::ostream & os, Sphere & sph)
		{
			return sph.printObject(os);
		}
	private:
		Eigen::Vector3d sphereCenter_;
		double sphereRadius_;
		std::string sphereColour_;
         	std::ostream & printObject(std::ostream & os); 
};

#endif // SPHERE_HPP
