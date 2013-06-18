#ifndef SPHERE_H
#define SPHERE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Config.h"

#include <Eigen/Dense>

/*!
 * \brief Class describing a sphere.
 * \author Roberto Di Remigio
 * \date 2011
 */

class Sphere
{
  	public:
		Sphere(){}
		Sphere(Eigen::Vector3d & _center, double _radius, const std::string & _colour = "Violet" )
			: sphereCenter(_center), sphereRadius(_radius), sphereColour(_colour) {}
		~Sphere(){}
                double getSphereRadius(){ return sphereRadius; }
                void setSphereRadius( double _radius ){ sphereRadius = _radius; }
		Eigen::Vector3d & getSphereCenter(){ return sphereCenter; }
                double getSphereCenter(int i){ return sphereCenter(i); }
                void setSphereCenter( Eigen::Vector3d & _coord ){ sphereCenter = _coord; }
		std::string & getSphereColour(){ return sphereColour; }
                void setSphereColour( std::string & _colour ){ sphereColour = _colour; }


         	friend inline void swap(Sphere & left, Sphere & right);
         	inline void swap(Sphere & other);
	        /// Assignment operator.
                Sphere & operator=(Sphere other);
                
                friend std::ostream& operator<<(std::ostream & o, Sphere & s);
	
	private:
		Eigen::Vector3d sphereCenter;
		double sphereRadius;
		std::string sphereColour;
         	std::ostream & printObject(std::ostream & os); 
 
};

#endif
