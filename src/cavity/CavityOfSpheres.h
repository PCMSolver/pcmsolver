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
         virtual VectorXd & getTessRadius(){return tessRadius;}
	 virtual int getNSpheres(){return nSpheres;}
     	 virtual void setNSpheres(int n){nSpheres = n;}
     	 virtual Matrix3Xd & getTessSphereCenter(){return tessSphereCenter;}
     	 virtual double getTessRadius(int i){return tessRadius(i);}
     	 virtual vector<Sphere> & getSpheres(){ return spheres; }
     	 virtual int getMode(){return mode;}
         virtual void setMode(const std::string & mode);
     	 virtual void setMode(int mode);

 private:
     	 enum SphereMode {Explicit, Atoms, Implicit};
	 SphereMode mode;
     	 int nSpheres;
     	 Matrix3Xd tessSphereCenter;
     	 VectorXd tessRadius;
     	 vector<Sphere> spheres;
};

#endif
