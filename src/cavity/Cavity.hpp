#ifndef CAVITY_HPP
#define CAVITY_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Sphere.hpp"

/*!
 *	\file Cavity.hpp
 *	\class Cavity
 *	\brief Abstract Base Class for cavities. 
 *	\author Krzysztof Mozgawa
 *	\date 2011 
 *
 * 	This class represents a cavity made of spheres, its surface being discretized in
 *      terms of finite elements.
 */

class Cavity
{
	public:
		//! Default constructor
		Cavity() : nElements(0), built(false) {}
		/*! \brief Constructor from spheres
		 *  \param[in] _spheres an STL vector containing the spheres making up the cavity.
		 */
		Cavity(const std::vector<Sphere> & _spheres) : spheres(_spheres), built(false)
			{
                		nSpheres = spheres.size();
				sphereCenter.resize(Eigen::NoChange, nSpheres);
				sphereRadius.resize(nSpheres);
				for (int i = 0; i < nSpheres; ++i) 
				{
					sphereCenter.col(i) = spheres[i].getSphereCenter();
					sphereRadius(i) = spheres[i].getSphereRadius();
				}
			}
                virtual ~Cavity(){}
		/*! \brief Creates the cavity and discretize its surface. 
		 */
                virtual void makeCavity() = 0;                                         
                virtual Eigen::Matrix3Xd & getElementCenter() { return elementCenter; }
                virtual Eigen::Vector3d getElementCenter(int i) { return elementCenter.col(i); }
                virtual Eigen::Matrix3Xd & getElementNormal() { return elementNormal; }
                virtual Eigen::Vector3d getElementNormal(int i) { return elementNormal.col(i); }
                virtual Eigen::VectorXd & getElementArea() { return elementArea; }
                virtual double getElementArea(int i) { return elementArea(i); }
                virtual int size() { return nElements; }
     	 	virtual std::vector<Sphere> & getSpheres() { return spheres; }
	  	virtual int getNSpheres() { return nSpheres; }
     	 	virtual void setNSpheres(int n) { nSpheres = n; }
                virtual Eigen::VectorXd & getSphereRadius() { return sphereRadius; }                
                virtual Eigen::Matrix3Xd & getSphereCenter() { return sphereCenter; }
                virtual Eigen::VectorXd & getElementRadius() { return elementRadius; }
     	        virtual double getElementRadius(int i) { return elementRadius(i); }
     	        virtual Eigen::Matrix3Xd & getElementSphereCenter() { return elementSphereCenter; }
	       	bool isBuilt() { return built; }
                double compPolarizationEnergy();
                double compPolarizationEnergy(const std::string & potential, const std::string & charge);
                friend std::ostream & operator<<(std::ostream & os, Cavity & cavity)
		{
			return cavity.printCavity(os);
		}
	 protected:
		std::vector<Sphere> spheres;
                int nElements;
                bool built;
		Eigen::Matrix3Xd elementCenter;
		Eigen::Matrix3Xd elementNormal;
		Eigen::VectorXd elementArea;
          	int nSpheres;
	 	Eigen::Matrix3Xd elementSphereCenter;
	 	Eigen::VectorXd elementRadius;
	        Eigen::Matrix3Xd sphereCenter;
	        Eigen::VectorXd sphereRadius;
	 private:
                virtual std::ostream & printCavity(std::ostream & os) = 0;  
};

#endif // CAVITY_HPP
