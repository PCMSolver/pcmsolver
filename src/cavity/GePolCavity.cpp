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

#include <boost/lexical_cast.hpp>

#include "Sphere.hpp"
#include "Symmetry.hpp"

extern "C" void generatecavity_cpp(int * maxts, int * maxsph, int * maxvert,
                             double * xtscor, double * ytscor, double * ztscor, double * ar, 
			     double * xsphcor, double * ysphcor, double * zsphcor, double * rsph, 
	                     int * nts, int * ntsirr, int * nesfp, int * addsph,
	                     double * xe, double * ye, double * ze, double * rin, 
			     double * avgArea, double * rsolv, double * ret, int * pgroup);

void GePolCavity::build(int maxts, int maxsph, int maxvert) 
{
       
	// This is a wrapper for the generatecavity_cpp function defined in the Fortran code PEDRA.
	// Here we allocate the necessary arrays to be passed to PEDRA, in particular we allow
	// for the insertion of additional spheres as in the most general formulation of the
	// GePol algorithm.

	double * xtscor  = new double[maxts];
	double * ytscor  = new double[maxts];
	double * ztscor  = new double[maxts];
	double * ar      = new double[maxts];
	double * xsphcor = new double[maxts];
	double * ysphcor = new double[maxts];
	double * zsphcor = new double[maxts];
	double * rsph    = new double[maxts];

	int nts = 0;
        int ntsirr = 0;

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
		
	double * xe = xv.data();
	double * ye = yv.data();
	double * ze = zv.data();

	double * rin = sphereRadius_.data();

        addedSpheres = 0;
	// Integer representing the point group
	int pg = pointGroup_.groupInteger();

        // Go PEDRA, Go!	
	generatecavity_cpp(&maxts, &maxsph, &maxvert,
                           xtscor, ytscor, ztscor, ar, xsphcor, ysphcor, zsphcor, rsph, 
		           &nts, &ntsirr, &nSpheres, &addedSpheres, 
			   xe, ye, ze, rin, &averageArea, &probeRadius, &minimalRadius, &pg);
	
	// The "intensive" part of updating the spheres related class data members will be of course
	// executed iff addedSpheres != 0
	if ( addedSpheres != 0 )
	{
		// Save the number of original spheres
 		int orig = nSpheres;
		// Update the nSpheres
		nSpheres += addedSpheres;
		// Resize sphereRadius and sphereCenter...
		sphereRadius.resize(nSpheres);
		sphereCenter.resize(Eigen::NoChange, nSpheres);
 		// Transfer radii and centers.
                // Eigen has no push_back function, so we need to traverse all the spheres...
		sphereRadius = sphereRadius_.head(nSpheres);
		for ( int i = 0; i < nSpheres; ++i )
		{
			sphereCenter(0, i) = xv(i);	
			sphereCenter(1, i) = yv(i);	
			sphereCenter(2, i) = zv(i);
		}
		// Now grow the vector<Sphere> containing the list of spheres
		for ( int i = orig;  i < nSpheres; ++i )
		{
			spheres.push_back(Sphere(sphereCenter.col(i), sphereRadius(i)));
		}
	}
    
        // Now take care of updating the rest of the cavity info.	
        nElements = static_cast<int>(nts);                                               
        nIrrElements = static_cast<int>(ntsirr);                                               
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
	// Check that no points are overlapping exactly
	// Do not perform float comparisons column by column. 
	// Instead form differences between columns and evaluate if they differ
	// from zero by more than a fixed threshold.
        // The indices of the equal elements are gathered in a std::pair and saved into a std::vector
	double threshold = 1.0e-12;
	std::vector< std::pair<int, int> > equal_elements;
        for( int i = 0; i < nElements; ++i )
	{
		for ( int j = i + 1; j < nElements; ++j)
		{
			Eigen::Vector3d difference = elementCenter.col(i) - elementCenter.col(j);
			if ( difference.isZero(threshold) )	
			{
				equal_elements.push_back(std::make_pair(i, j));
			}
		}
	}
	if ( equal_elements.size() != 0 )
	{
		// Not sure that printing the list of pairs is actually of any help...
		std::string list_of_pairs;
		for ( size_t i = 0; i < equal_elements.size(); ++i)
		{
			list_of_pairs += "(" + boost::lexical_cast<std::string>(equal_elements[i].first) 
				     + ", " + boost::lexical_cast<std::string>(equal_elements[i].second) + ")\n";
		}
		// Prepare the error message:
		std::string message = boost::lexical_cast<std::string>(equal_elements.size()) + " cavity finite element centers overlap exactly!\n" + list_of_pairs; 
		throw std::runtime_error(message);
	}
        // Calculate normal vectors 
        elementNormal = elementCenter - elementSphereCenter;
        for( int i = 0; i < nElements; ++i)
	{
    		elementNormal.col(i) /= elementNormal.col(i).norm();
    	}
   
        // Clean-up 
	delete[] xtscor;
	delete[] ytscor;
	delete[] ztscor;
	delete[] ar;
	delete[] xsphcor;
	delete[] ysphcor;
	delete[] zsphcor;
	delete[] rsph;
	
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
	os << "Probe radius = " << probeRadius << std::endl;
        if ( addedSpheres != 0 )
	{
		os << "Addition of extra spheres enabled" << std::endl;
 	}
	os << "Number of spheres = " << nSpheres << " [initial = " << nSpheres - addedSpheres << "; added = " << addedSpheres << "]" << std::endl;
        os << "Number of finite elements = " << nElements << std::endl;
	if (pointGroup_.groupInteger() != 0)
	{
		os << "Number of irreducible finite elements = " << nIrrElements;
	}
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
