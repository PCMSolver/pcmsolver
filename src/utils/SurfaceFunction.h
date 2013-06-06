#ifndef SURFACE_FUNCTION
#define SURFACE_FUNCTION

#include <iostream>
#include <string>
#include <map>

#include "Config.h"

#include <Eigen/Dense>

/** 

A basic surface function class
written by L. Frediani 2012

*/

class SurfaceFunction
{
 typedef std::map<std::string, SurfaceFunction *> SurfaceFunctionMap;

 public:
         SurfaceFunction() : nPoints(0), allocated(false) {}
         SurfaceFunction(const std::string & name_) : name(name_), allocated(false) {}
         SurfaceFunction(const std::string & name_, int nPoints_) : name(name_), nPoints(nPoints_) 
	 {
			values = Eigen::VectorXd::Zero(nPoints);
			allocated = true;
			// The SurfaceFunction registers itself in the SurfaceFunctionMap
//			Register();
	 } 							       
         SurfaceFunction(const std::string & name_, int nPoints_, double * values_) : name(name_), nPoints(nPoints_)            	
	 {
			values = Eigen::VectorXd::Zero(nPoints);
	 		allocated = true;
	 		for (int i = 0; i < nPoints; ++i)
	 		{
	 			values(i) = values_[i];
	 		}
	 		// The SurfaceFunction registers itself in the SurfaceFunctionMap
//	 		Register();
	 }
         ~SurfaceFunction()                                                                                                     
         {
	     std::cout << "Calling DTOR" << std::endl;
             allocated = false;
             // Upon destruction the SurfaceFunction unregisters itself from the SurfaceFunctionMap.
             // This should be enough to avoid dangling pointers...
//	     unRegister();
         }
                                                                                                                         
         /// Copy constructor
         SurfaceFunction(const SurfaceFunction & other) : name(other.name), nPoints(other.nPoints), values(other.values)
         {
             allocated = true;
             // The SurfaceFunction registers itself in the SurfaceFunctionMap
             // static member of the Cavity class.
  //           Register();
         }
                                                                                                                        
	 static SurfaceFunctionMap & initSurfaceFunctionMap()
	 {
	 	static SurfaceFunctionMap func;
	 	return func;
	 }	
                                                                                                                                 
         friend inline void swap(SurfaceFunction & left, SurfaceFunction & right);
         inline void swap(SurfaceFunction & other);
         /// Assignment operator
         SurfaceFunction & operator=(SurfaceFunction other);
         /// Multiplication operator: product of two SurfaceFunctions version (scalar product of the values vectors)
         double operator*(const SurfaceFunction & other);
         /// Addition-assignment operator
         SurfaceFunction & operator+=(const SurfaceFunction & other);
         /// Subtraction-assignment operator
         SurfaceFunction & operator-=(const SurfaceFunction & other);
         /// Multiplication-assignment operator. Defined only for the uniform scaling case of operator*
         SurfaceFunction & operator*=(double scaling);
                                                                                                                         
         std::string & getName(){ return name; }
         int getNPoints(){ return nPoints; }
         void setValue(int index_, double value_) { values(index_) = value_; }
         double getValue(int index_) {return values(index_);}
         Eigen::VectorXd & getVector(){ return values; }
         void allocate(int nPoints_){ values.resize(nPoints_); }
         bool isAllocated() { return allocated; }
         void clear();
                                                                                                                         
         void setValues(double * value);
         void getValues(double * value);
                                                                                                                         
         bool Register();
         bool unRegister();
                                                                                                                         
         friend std::ostream & operator<<(std::ostream & o, SurfaceFunction & s);

 private:
         virtual std::ostream & printObject(std::ostream & os); 
         std::string name;
         int nPoints;
         Eigen::VectorXd values;
         bool allocated;
         bool registered;
	 static SurfaceFunctionMap functions;
};


/// Addition operator
inline SurfaceFunction operator+(SurfaceFunction left, const SurfaceFunction & right)
{
	left += right;
	return left;
}

/// Subtraction operator
inline SurfaceFunction operator-(SurfaceFunction left, const SurfaceFunction & right)
{
	left -= right;
	return left;
}

/// Multiplication operator: uniform scaling of SurfaceFunction version
inline SurfaceFunction operator*(double scaling, SurfaceFunction & object)
{
	object *= scaling;
	return object;
}

#endif
