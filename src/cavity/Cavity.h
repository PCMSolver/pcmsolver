#ifndef CAVITY
#define CAVITY

#include <iostream>
#include <string>
#include <map>

#include "Config.h"

#include <Eigen/Dense>

#include "Sphere.h"

/** 

A basic cavity class
written by Krzysztof Mozgawa, 2011


*/

class Cavity
{
	public:
		Cavity() : nElements(0), built(false) {}
		Cavity(const std::vector<Sphere> & _spheres) : spheres(_spheres)
			{
				spheres = _spheres;
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
                virtual void makeCavity() = 0;                                         
                virtual void writeOutput(std::string &filename);
		
		/// Functions related to the finite elements making up the cavity
                virtual Eigen::Matrix3Xd & getElementCenter(){return elementCenter;}
                virtual Eigen::Vector3d getElementCenter(int i){return elementCenter.col(i);}
                virtual Eigen::Matrix3Xd & getElementNormal(){return elementNormal;}
                virtual Eigen::Vector3d getElementNormal(int i){return elementNormal.col(i);}
                virtual Eigen::VectorXd & getElementArea(){return elementArea;}
                virtual double getElementArea(int i){return elementArea(i);}
                virtual int size(){return nElements;}
		
		/// Functions related to the spheres composing the cavity surface
     	 	virtual vector<Sphere> & getSpheres(){ return spheres; }
	  	virtual int getNSpheres(){return nSpheres;}
     	 	virtual void setNSpheres(int n){nSpheres = n;}
                virtual Eigen::VectorXd & getSphereRadius(){return sphereRadius;}                
                virtual Eigen::Matrix3Xd & getSphereCenter(){return sphereCenter;}
                virtual Eigen::VectorXd & getElementRadius(){return elementRadius;}
     	        virtual double getElementRadius(int i){return elementRadius(i);}
     	        virtual Eigen::Matrix3Xd & getElementSphereCenter(){return elementSphereCenter;}
     	        virtual int getMode(){return mode;}
                virtual void setMode(const std::string & mode);
     	        virtual void setMode(int mode);
               
	       	bool isBuilt(){return built;}
                double compPolarizationEnergy();
                double compPolarizationEnergy(const std::string & potential, const std::string & charge);
               
     	        enum SphereMode {Explicit, Atoms, Implicit};
                enum chargeType{Nuclear, Electronic};
                                                                                       
                friend std::ostream& operator<<(std::ostream & o, Cavity & c);
    
	 protected:
                virtual std::ostream & printObject(std::ostream & os);  
                int nElements;
                bool built;
		Eigen::Matrix3Xd elementCenter;
		Eigen::Matrix3Xd elementNormal;
		Eigen::VectorXd elementArea;
	        SphereMode mode;
          	int nSpheres;
	 	Eigen::Matrix3Xd elementSphereCenter;
	 	Eigen::VectorXd elementRadius;
	        Eigen::Matrix3Xd sphereCenter;
	        Eigen::VectorXd sphereRadius;
     	        vector<Sphere> spheres;
};


#endif
