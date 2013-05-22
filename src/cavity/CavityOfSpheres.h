#ifndef CAVITYOFSPHERES_H
#define CAVITYOFSPHERES_H 

#include <iostream>
#include <string>
#include <vector>

#include <Config.h>

#include <Eigen/Dense>

#include "Cavity.h"
#include "Sphere.h"


class CavityOfSpheres : public Cavity 
{
 public:
	 CavityOfSpheres(){}
	 ~CavityOfSpheres(){}
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
     	 enum SphereMode {Explicit, Atoms, Implicit};

 protected:
	 SphereMode mode;
     	 int nSpheres;
	 Eigen::Matrix3Xd elementSphereCenter;
	 Eigen::VectorXd elementRadius;
	 Eigen::Matrix3Xd sphereCenter;
	 Eigen::VectorXd sphereRadius;
     	 vector<Sphere> spheres;
};

#endif
