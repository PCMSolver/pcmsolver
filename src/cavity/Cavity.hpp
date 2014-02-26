#ifndef CAVITY_HPP
#define CAVITY_HPP

#include <iosfwd>
#include <vector>

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

#include "Sphere.hpp"
#include "Symmetry.hpp"

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
	protected:
		std::vector<Sphere> spheres;
                int nElements;
		int nIrrElements;
                bool built;
		Eigen::Matrix3Xd elementCenter;
		Eigen::Matrix3Xd elementNormal;
		Eigen::VectorXd elementArea;
          	int nSpheres;
	 	Eigen::Matrix3Xd elementSphereCenter;
	 	Eigen::VectorXd elementRadius;
	        Eigen::Matrix3Xd sphereCenter;
	        Eigen::VectorXd sphereRadius;
		Symmetry pointGroup_;
	private:
		/*! \brief Creates the cavity and discretizes its surface. 
		 */
                virtual void makeCavity() = 0; 
                virtual std::ostream & printCavity(std::ostream & os) = 0;  
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
					sphereCenter.col(i) = spheres[i].sphereCenter();
					sphereRadius(i) = spheres[i].sphereRadius();
				}
			}
		Cavity(const std::vector<Sphere> & _spheres, int pg) : spheres(_spheres), built(false), pointGroup_(buildGroup(pg))
			{
                		nSpheres = spheres.size();
				sphereCenter.resize(Eigen::NoChange, nSpheres);
				sphereRadius.resize(nSpheres);
				for (int i = 0; i < nSpheres; ++i) 
				{
					sphereCenter.col(i) = spheres[i].sphereCenter();
					sphereRadius(i) = spheres[i].sphereRadius();
				}
			}
                virtual ~Cavity() {}
                Eigen::Matrix3Xd & getElementCenter() { return elementCenter; }
                Eigen::Vector3d getElementCenter(int i) { return elementCenter.col(i); }
                Eigen::Matrix3Xd & getElementNormal() { return elementNormal; }
                Eigen::Vector3d getElementNormal(int i) { return elementNormal.col(i); }
                Eigen::VectorXd & getElementArea() { return elementArea; }
                double getElementArea(int i) { return elementArea(i); }
                int size() { return nElements; }
		int irreducible_size() { return nIrrElements; }
		virtual Symmetry pointGroup() const { return pointGroup_; } 
     	 	std::vector<Sphere> & getSpheres() { return spheres; }
	  	int getNSpheres() { return nSpheres; }
     	 	void setNSpheres(int n) { nSpheres = n; }
                Eigen::VectorXd & getSphereRadius() { return sphereRadius; }                
                Eigen::Matrix3Xd & getSphereCenter() { return sphereCenter; }
                Eigen::VectorXd & getElementRadius() { return elementRadius; }
     	        double getElementRadius(int i) { return elementRadius(i); }
     	        Eigen::Matrix3Xd & getElementSphereCenter() { return elementSphereCenter; }
		/*! \brief Save cavity specification to file.
		 *
		 *  The cavity specification contains:
		 *   0. the number of finite elements, nElements;
		 *   1. the weight of the finite elements, elementArea;
		 *   2. the radius of the finite elements, elementRadius;
		 *   3. the centers of the finite elements, elementCenter;
		 *   4. the normal vectors relative to the centers, elementNormal.
		 *  Each of these objects is saved in a separate .npy binary file
		 *  and compressed into one .npz file.
		 *  Notice that this is just the minimal set of data needed to 
		 *  restart an energy calculation.
		 */
		virtual void saveCavity(const std::string & fname = "cavity.npz");
		/*! \brief Load cavity specification from file.
		 */
		virtual void loadCavity(const std::string & fname = "cavity.npz");
	       	bool isBuilt() { return built; }
                friend std::ostream & operator<<(std::ostream & os, Cavity & cavity)
		{
			return cavity.printCavity(os);
		}
};

#endif // CAVITY_HPP
