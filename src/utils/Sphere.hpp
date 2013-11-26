#ifndef SPHERE_HPP
#define SPHERE_HPP

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

/*! \file Sphere.hpp
 *  \class Sphere
 *  \brief Class describing a sphere.
 *  \author Roberto Di Remigio
 *  \date 2011
 */

class Sphere
{
  	public:
		Sphere() {}
		Sphere(Eigen::Vector3d & _center, double _radius, const std::string & _colour = "Violet" )
			: sphereCenter(_center), sphereRadius(_radius), sphereColour(_colour) {}
		~Sphere() {}
                double getSphereRadius() const { return sphereRadius; }
                void setSphereRadius( double _radius ) { sphereRadius = _radius; }
		Eigen::Vector3d & getSphereCenter() { return sphereCenter; }
                double getSphereCenter(int i) const { return sphereCenter(i); }
                void setSphereCenter( Eigen::Vector3d & _coord ){ sphereCenter = _coord; }
		std::string & getSphereColour() { return sphereColour; }
                void setSphereColour( std::string & _colour ){ sphereColour = _colour; }
         	friend inline void swap(Sphere & left, Sphere & right);
         	inline void swap(Sphere & other);
	        /// Assignment operator.
                Sphere & operator=(Sphere other);
                friend std::ostream& operator<<(std::ostream & os, Sphere & sph)
		{
			return sph.printObject(os);
		}
	private:
		Eigen::Vector3d sphereCenter;
		double sphereRadius;
		std::string sphereColour;
         	std::ostream & printObject(std::ostream & os); 
};

#endif // SPHERE_HPP
