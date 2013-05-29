#ifndef CAVITY
#define CAVITY

#include <iostream>
#include <string>
#include <map>

#include "Config.h"

#include <Eigen/Dense>


/** 

A basic cavity class
written by Krzysztof Mozgawa, 2011


*/

class SurfaceFunction;


class Cavity
{
	protected:
		typedef std::pair< std::string, SurfaceFunction * > SurfaceFunctionPair;
		typedef std::map< std::string, SurfaceFunction * > SurfaceFunctionMap;
	public:
		Cavity() : nElements(0), built(false) {}
                virtual ~Cavity(){}
                virtual void makeCavity() = 0;                                         
                virtual void writeOutput(std::string &filename);
                virtual Eigen::Matrix3Xd & getElementCenter(){return elementCenter;}
                virtual Eigen::Vector3d getElementCenter(int i){return elementCenter.col(i);}
                virtual Eigen::Matrix3Xd & getElementNormal(){return elementNormal;}
                virtual Eigen::Vector3d getElementNormal(int i){return elementNormal.col(i);}
                virtual Eigen::VectorXd & getElementArea(){return elementArea;}
                virtual double getElementArea(int i){return elementArea(i);}
                virtual int size(){return nElements;}
                bool isBuilt(){return built;}
                double compPolarizationEnergy();
                double compPolarizationEnergy(const std::string & potential, const std::string & charge);
                void appendNewFunction(const std::string & name);
                void setFunction(const std::string & name, double * values);
                SurfaceFunction & getFunction(const std::string & name);
                bool functionExists(const std::string & name) 
		{
			SurfaceFunctionMap::const_iterator i = functions.find(name);
			return i != functions.end();	
                }
                enum chargeType{Nuclear, Electronic};
                                                                                       
                friend std::ostream& operator<<(std::ostream & o, Cavity & c);
    
	 protected:
                virtual std::ostream & printObject(std::ostream & os);  
                int nElements;
                bool built;
		Eigen::Matrix3Xd elementCenter;
		Eigen::Matrix3Xd elementNormal;
		Eigen::VectorXd elementArea;
                SurfaceFunctionMap functions;
};


#endif
