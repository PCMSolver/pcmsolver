#include "GePolCavity.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Config.hpp"
#include "FCMangle.hpp"

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

extern "C" void generatecavity_cpp(double *xtscor, double *ytscor, double *ztscor, double *ar, double *xsphcor, double *ysphcor, 
			     double *zsphcor, double *rsph, int *nts, int *nesfp, double *xe, double *ye, double *ze, double *rin, 
			     double *avgArea, double *rsolv, double * ret, int * pgroup, double* work, int* lwork);

void GePolCavity::makeCavity()
{
	makeCavity(10000, 10000000);
}

void GePolCavity::makeCavity(int maxts, int lwork) 
{
       
	// This is a wrapper for the generatecavity_cpp function defined in the Fortran code PEDRA.
	// Here we allocate the necessary arrays to be passed to PEDRA, in particular we allow
	// for the insertion of additional spheres as in the most general formulation of the
	// GePol algorithm.

	double *xtscor  = new double[maxts];
	double *ytscor  = new double[maxts];
	double *ztscor  = new double[maxts];
	double *ar      = new double[maxts];
	double *xsphcor = new double[maxts];
	double *ysphcor = new double[maxts];
	double *zsphcor = new double[maxts];
	double *rsph    = new double[maxts];
	double *work    = new double[lwork];

	int nts = 0;

        // If there's an overflow in the number of spheres PEDRA will die.
        // The maximum number of spheres in PEDRA is set to 200 (primitive+additional)
        // so the integer here declared is just to have enough space C++ side to pass everything back.
	int maxAddedSpheres = 200;
	
	// Allocate vectors of size equal to nSpheres + maxAddedSpheres where maxAddedSpheres is the 
	// maximum number of spheres we allow the algorithm to add to our original set.
	// If this number is exceeded, then the algorithm crashes (should look into this...)
	// After the cavity is generated we will update ALL the class data members, both related
	// to spheres and finite elements so that the cavity is fully formed.
	
	Eigen::VectorXd xv = Eigen::VectorXd::Zero(nSpheres + maxAddedSpheres);
	Eigen::VectorXd yv = Eigen::VectorXd::Zero(nSpheres + maxAddedSpheres);
	Eigen::VectorXd zv = Eigen::VectorXd::Zero(nSpheres + maxAddedSpheres);
	Eigen::VectorXd sphereRadius_ = Eigen::VectorXd::Zero(nSpheres + maxAddedSpheres); // Not to be confused with the data member inherited from Cavity!!!
	
	for ( int i = 0; i < nSpheres; ++i ) 
	{
		for ( int j = 0; j < 3; ++j )
		{
			xv(i) = sphereCenter(0, i);
			yv(i) = sphereCenter(1, i);
			zv(i) = sphereCenter(2, i);
		}
		sphereRadius_(i) = sphereRadius(i);
	}
		
	double *xe = xv.data();
	double *ye = yv.data();
	double *ze = zv.data();

	double *rin = sphereRadius_.data();

	int pointGroup = 0;

        // Go PEDRA, Go!	
	generatecavity_cpp(xtscor, ytscor, ztscor, ar, xsphcor, ysphcor, zsphcor, rsph, &nts, &nSpheres, 
						xe, ye, ze, rin, &averageArea, &probeRadius, &minimalRadius, &pointGroup, work, &lwork);
	
	// We now create come Eigen temporaries to be used in the post-processing of the spheres data
	// coming out from generatecavity_cpp_
	
	Eigen::VectorXd rtmp = Eigen::VectorXd::Zero(nSpheres + maxAddedSpheres);
	Eigen::VectorXd xtmp = Eigen::VectorXd::Zero(nSpheres + maxAddedSpheres);
	Eigen::VectorXd ytmp = Eigen::VectorXd::Zero(nSpheres + maxAddedSpheres);
	Eigen::VectorXd ztmp = Eigen::VectorXd::Zero(nSpheres + maxAddedSpheres);

	// The first nSpheres elements of these temporaries will still be those of the original set of spheres	
	rtmp = sphereRadius_.head(nSpheres);
	xtmp = xv.head(nSpheres);
	ytmp = yv.head(nSpheres);
	ztmp = zv.head(nSpheres);
	// Traverse the sphereRadius vector (starting from the nSpheres index) and count the number of 
	// additional spheres that the algorithm created (characterized by non-zero radius)
	addedSpheres = 0;
	for ( int i = nSpheres; i < nSpheres + maxAddedSpheres; i++ ) 
	{
		if ( sphereRadius_(i) != 0.0) 
		{
			rtmp(i) = sphereRadius_(i);
			xtmp(i) = xv(i);
			ytmp(i) = yv(i);
			ztmp(i) = zv(i);
			addedSpheres += 1;
		}
	}

	// The "intensive" part of updating the spheres related class data members will be of course
	// executed iff addedSpheres != 0
	if ( addedSpheres != 0 )
	{
		std::cout << "The PEDRA algorithm added " << addedSpheres << " new spheres to the original list." << std::endl;
		// First of all update the nSpheres
		nSpheres += addedSpheres;
		// Resize sphereRadius and sphereCenter...
		sphereRadius.resize(nSpheres);
		sphereCenter.resize(Eigen::NoChange, nSpheres);
		// ...clear vector<Sphere> spheres...
		spheres.clear();
		// ...and update their content
		for ( int i = 0; i < nSpheres; ++i )
		{
			sphereRadius(i) = rtmp(i);
			for ( int j = 0; j < 3; ++j )
			{
				sphereCenter(0, i) = xtmp(i);
				sphereCenter(1, i) = ytmp(i);
				sphereCenter(2, i) = ztmp(i);
			}
			Eigen::Vector3d cent = sphereCenter.col(i);
			Sphere sph(cent, sphereRadius(i));
			spheres.push_back(sph);
		}
	}
    	
        nElements = int(nts);                                               
        elementCenter.resize(Eigen::NoChange, nElements);
        elementSphereCenter.resize(Eigen::NoChange, nElements);
        elementNormal.resize(Eigen::NoChange, nElements);
        elementArea.resize(nElements);
        elementRadius.resize(nElements);
        for( int i = 0; i < nElements; ++i )
	{
    		elementCenter(0,i) = xtscor[i];
    		elementCenter(1,i) = ytscor[i];
    		elementCenter(2,i) = ztscor[i];
    		elementArea(i) = ar[i];
    		elementSphereCenter(0,i) = xsphcor[i];
    		elementSphereCenter(1,i) = ysphcor[i];
    		elementSphereCenter(2,i) = zsphcor[i];
    		elementRadius(i) = rsph[i];
        }
    
        elementNormal = elementCenter - elementSphereCenter;
        for( int i = 0; i < nElements; ++i)
	{
    		elementNormal.col(i) /= elementNormal.col(i).norm();
    	}
    
	delete[] xtscor;
	delete[] ytscor;
	delete[] ztscor;
	delete[] ar;
	delete[] xsphcor;
	delete[] ysphcor;
	delete[] zsphcor;
	delete[] rsph;
	delete[] work;
	
	built = true;
}

std::ostream & GePolCavity::printCavity(std::ostream & os) 
{
	/*
	  We should print the cavity.off file here, just to
	  get a prettier cavity image.
	*/
        os << "Cavity type: GePol" << std::endl;
	os << "Average area = " << averageArea << " AU^2" << std::endl;
	os << "Addition of extra spheres enabled" << std::endl;
	os << "Probe radius = " << probeRadius << std::endl;
	os << "Number of spheres = " << nSpheres << " [initial = " << nSpheres - addedSpheres << "; added = " << addedSpheres << "]" << std::endl;
        os << "Number of finite elements = " << nElements;
        /*for(int i = 0; i < nElements; i++) 
	{
		os << std::endl;
		os << i+1 << " ";
		os << elementCenter(0,i) << " ";
		os << elementCenter(1,i) << " ";
		os << elementCenter(2,i) << " ";
		os << elementArea(i) << " ";
        }*/
	return os;
}
