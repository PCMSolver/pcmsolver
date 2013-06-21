#ifndef GREENSFUNCTIONSUM_H
#define GREENSFUNCTIONSUM_H

#include "Config.h"

/*! \file GreensFunctionSum.h
 *  \class GreensFunctionSum
 *  \brief A class to describe the sum of two Green's functions
 *  \author Luca Frediani
 *  \date 2012
 *
 * I am for now limiting myself to one template parameter. 
 * In principle three templates are necessary here. 
 * The question is whether it is worth implementng the full version.
 * RDR couldn't we implement overloaded versions of operator+ and operator-
 * instead of having an ad hoc class?
 */


template<class T>
class GreensFunctionSum : public GreensFunction<T>
{
	public:
		GreensFunctionSum(){}
		GreensFunctionSum(GreensFunction<T> &first, GreensFunction<T> &second);
//              GreensFunctionSum(Section green);
                ~GreensFunctionSum(){delete greenFirst; delete greenSecond;}   
                double evalf(Eigen::Vector3d &p1, Eigen::Vector3d &p2);
                double evald(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2);
                double compDiagonalElementS(double area);
                double compDiagonalElementD(double area, double radius);
	
	protected:
                T evalGreensFunction(T * source, T * probe);
                GreensFunctionInterface* greenFirst;
                GreensFunctionInterface* greenSecond;
};

#endif // GREENSFUNCTIONSUM_H
