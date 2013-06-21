#ifndef METALSPHERE_H
#define METALSPHERE_H

#include <complex>

#include "Config.h"

/*! \file MetalSphere.h
 *  \class MetalSphere
 *  \brief Class to describe spherical metal nanoparticles.
 *  \author Luca Frediani
 *  \date 2011
 *
 *  This class is a wrapper around the FORTRAN routines written
 *  by Stefano Corni et al. to take into account the presence
 *  of metal nanoparticles.
 */


class MetalSphere : public GreensFunction<double>
{
	private:
		typedef std::complex<double> dcomplex;

	public:
		MetalSphere(double eps, double epsRe, double epsIm, Eigen::Vector3d &pos, double radius)
			: epsSolvent(eps), epsMetal(dcomplex(epsRe, epsIm)), sphPosition(pos), sphRadius(radius) 
		{
			uniformFlag = false;
		}
//    MetalSphere(Section green);
                ~MetalSphere(){};                                              
                double evald(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2);
                double compDiagonalElementS(double area);
                double compDiagonalElementD(double area, double radius);
 	
	private:
                double evalGreensFunction(double * source, double * probe); 
                double epsSolvent;
                dcomplex epsMetal;
                Eigen::Vector3d sphPosition;
                double sphRadius;
};
#endif // METALSPHERE_H
