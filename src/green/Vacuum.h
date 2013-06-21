#ifndef VACUUM_H
#define VACUUM_H

#include <Eigen/Dense>

#include "Config.h"

#include "GreensFunction.h"

/*! \file Vacuum.h
 *  \class Vacuum
 *  \brief Class for the Green's function of vacuum.
 *  \author Luca Frediani
 *  \date 2011
 *
 *  This class represents the Green's function for the vacuum. 
 *
 */

template<typename T>
class Vacuum : public GreensFunction<T>
{
	public:
		Vacuum(){GreensFunction<T>::uniformFlag = true;};
    		~Vacuum(){};
    		double evald(Eigen::Vector3d & direction, Eigen::Vector3d & p1, Eigen::Vector3d & p2);
    		double compDiagonalElementS(double area);
    		double compDiagonalElementD(double area, double radius);
    		double getDielectricConstant(){return 1.0;}

 	private:
    		T evalGreensFunction(T * source, T * probe);
};

namespace
{
}

#endif
