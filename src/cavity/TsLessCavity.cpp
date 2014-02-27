#include "TsLessCavity.hpp"

#include <iostream>
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

extern "C" void generate_tslesscavity_cpp(double *xtscor, double *ytscor, double *ztscor, double *ar, double *xsphcor, double *ysphcor, 
			     double *zsphcor, double *rsph, int *nts, int *nesfp, double *xe, double *ye, double *ze, double *rin, 
			     double *avgArea, double *rsolv, double * dmin, int * nord, double* work, int* lwork);

void TsLessCavity::build(int maxts, int maxsph, int maxvert) 
{
	// This is a wrapper for the generatecavity_cpp_ function defined in the Fortran code PEDRA.
	// Here we allocate the necessary arrays to be passed to PEDRA, in particular we allow
	// for the insertion of additional spheres as in the most general formulation of the
	// GePol algorithm.

	int lwork = 1000000;
	double *xtscor  = new double[maxts];
	double *ytscor  = new double[maxts];
	double *ztscor  = new double[maxts];
	double *ar      = new double[maxts];
	double *xsphcor = new double[maxts];
	double *ysphcor = new double[maxts];
	double *zsphcor = new double[maxts];
	double *rsph    = new double[maxts];
	double *work    = new double[lwork];
	
	int nts;
	int maxAddedSpheres = 200;
	
	// Allocate vectors of size equal to nSpheres_ + maxAddedSpheres where maxAddedSpheres is the 
	// maximum number of spheres we allow the algorithm to add to our original set.
	// If this number is exceeded, then the algorithm crashes (should look into this...)
	// After the cavity is generated we will update ALL the class data members, both related
	// to spheres and finite elements so that the cavity is fully formed.
	
	Eigen::VectorXd xv = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
	Eigen::VectorXd yv = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
	Eigen::VectorXd zv = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
	Eigen::VectorXd radii_scratch = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres); // Not to be confused with the data member inherited from Cavity!!!
	
	for ( int i = 0; i < nSpheres_; ++i ) 
	{
		for ( int j = 0; j < 3; ++j )
		{
			xv(i) = sphereCenter_(0, i);
			yv(i) = sphereCenter_(1, i);
			zv(i) = sphereCenter_(2, i);
		}
		radii_scratch(i) = sphereRadius_(i);
	}
		
	double *xe = xv.data();
	double *ye = yv.data();
	double *ze = zv.data();

	double *rin = radii_scratch.data();

        // Go TsLess, Go!	
	generate_tslesscavity_cpp(xtscor, ytscor, ztscor, ar, xsphcor, ysphcor, zsphcor, rsph, &nts, &nSpheres_, 
						xe, ye, ze, rin, &averageArea, &probeRadius, &minDistance, &derOrder, work, &lwork);
        throw std::runtime_error("TsLessCavity Fortran backend not yet implemented!");
	
	// We now create come Eigen temporaries to be used in the post-processing of the spheres data
	// coming out from generatecavity_cpp_
	
	Eigen::VectorXd rtmp = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
	Eigen::VectorXd xtmp = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
	Eigen::VectorXd ytmp = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
	Eigen::VectorXd ztmp = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);

	// The first nSpheres_ elements of these temporaries will still be those of the original set of spheres	
	rtmp = radii_scratch.head(nSpheres_);
	xtmp = xv.head(nSpheres_);
	ytmp = yv.head(nSpheres_);
	ztmp = zv.head(nSpheres_);
	// Traverse the sphereRadius vector (starting from the nSpheres_ index) and count the number of 
	// additional spheres that the algorithm created (characterized by non-zero radius)
	addedSpheres = 0;
	for ( int i = nSpheres_; i < nSpheres_ + maxAddedSpheres; i++ ) 
	{
		if ( radii_scratch(i) != 0.0) 
		{
			rtmp(i) = radii_scratch(i);
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
		// First of all update the nSpheres_
		nSpheres_ += addedSpheres;
		// Resize sphereRadius and sphereCenter...
		sphereRadius_.resize(nSpheres_);
		sphereCenter_.resize(Eigen::NoChange, nSpheres_);
		// ...clear vector<Sphere> spheres...
		spheres_.clear();
		// ...and update their content
		for ( int i = 0; i < nSpheres_; ++i )
		{
			sphereRadius_(i) = rtmp(i);
			for ( int j = 0; j < 3; ++j )
			{
				sphereCenter_(0, i) = xtmp(i);
				sphereCenter_(1, i) = ytmp(i);
				sphereCenter_(2, i) = ztmp(i);
			}
			Eigen::Vector3d cent = sphereCenter_.col(i);
			Sphere sph(cent, sphereRadius_(i));
			spheres_.push_back(sph);
		}
	}
    	
        nElements_ = int(nts);                                               
        elementCenter_.resize(Eigen::NoChange, nElements_);
        elementSphereCenter_.resize(Eigen::NoChange, nElements_);
        elementNormal_.resize(Eigen::NoChange, nElements_);
        elementArea_.resize(nElements_);
        elementRadius_.resize(nElements_);
        for( int i = 0; i < nElements_; ++i )
	{
    		elementCenter_(0,i) = xtscor[i];
    		elementCenter_(1,i) = ytscor[i];
    		elementCenter_(2,i) = ztscor[i];
    		elementArea_(i) = ar[i];
    		elementSphereCenter_(0,i) = xsphcor[i];
    		elementSphereCenter_(1,i) = ysphcor[i];
    		elementSphereCenter_(2,i) = zsphcor[i];
    		elementRadius_(i) = rsph[i];
        }
    
        elementNormal_ = elementCenter_ - elementSphereCenter_;
        for( int i = 0; i < nElements_; ++i)
	{
    		elementNormal_.col(i) /= elementNormal_.col(i).norm();
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

std::ostream & TsLessCavity::printCavity(std::ostream & os) 
{
	/*
	  We should print the cavity.off file here, just to
	  get a prettier cavity image.
	*/
        os << "Cavity type: GePol" << std::endl;
	os << "Average area = " << averageArea << " AU^2" << std::endl;
	os << "Addition of extra spheres enabled" << std::endl;
	os << "Probe radius = " << probeRadius << std::endl;
	os << "Number of spheres = " << nSpheres_ << " [initial = " << nSpheres_ - addedSpheres << "; added = " << addedSpheres << "]" << std::endl;
        os << "Number of finite elements = " << nElements_;
        /*for(int i = 0; i < nElements_; i++) 
	{
		os << std::endl;
		os << i+1 << " ";
		os << elementCenter_(0,i) << " ";
		os << elementCenter_(1,i) << " ";
		os << elementCenter_(2,i) << " ";
		os << elementArea_(i) << " ";
        }*/
	return os;
}
