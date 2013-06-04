/*

  GePol c++ interface and wrapper methods
  written by Krzysztof Mozgawa, 2011

*/
#include <iostream>
#include <fstream> 
#include <string>
#include <vector>

#include <Eigen/Dense>

using std::cout;
using std::endl;

#include "GePolCavity.h"
#include "CavityFactory.h"

bool GePolCavity::registered = false;

bool GePolCavity::Register()
{
	std::cout << "Register GePolCavity" << std::endl;
	GePolCavity::registered = CavityFactory::TheCavityFactory().registerCavity("GePol", GePolCavity::Create);
	return GePolCavity::registered;
}

void GePolCavity::writeOutput(string &filename){

    ofstream output;
    output.open(filename.c_str(), fstream::out);

    output << nElements << endl;
    for(int i=0; i < nElements; i++) {
		output << elementCenter(0,i) << " ";
		output << elementCenter(1,i) << " ";
		output << elementCenter(2,i) << " ";
		output << elementArea(i) << " ";
		output << elementSphereCenter(0,i) << " ";
		output << elementSphereCenter(1,i) << " ";
		output << elementSphereCenter(2,i) << " ";
		output << elementRadius(i) << endl;
    }
    output.close();
}

extern "C" {
    void generatecavity_cpp_(double *xtscor, double *ytscor, double *ztscor, 
							 double *ar, double *xsphcor, double *ysphcor, 
							 double *zsphcor, double *rsph, int *nts, int *nesfp,
							 double *xe, double *ye, double *ze, double *rin, 
							 double *avgArea, double *rsolv, double* work, int* lwork);
}

void GePolCavity::makeCavity(){
	makeCavity(10000, 10000000);
}

void GePolCavity::makeCavity(int maxts, int lwork) {
       
	// This is a wrapper for the generatecavity_cpp_ function defined in the Fortran code PEDRA.
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
	
	int nts;
	
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

        // Go PEDRA, Go!	
	generatecavity_cpp_(xtscor, ytscor, ztscor, ar, xsphcor, ysphcor, zsphcor, rsph, &nts, &nSpheres, 
						xe, ye, ze, rin, &averageArea, &probeRadius, work, &lwork);
	
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
			Vector3d cent = sphereCenter.col(i);
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

ostream & GePolCavity::printObject(ostream & os) {
	/*
	  We should print the cavity.off file here, just to
	  get a prettier cavity image.
	*/
	os << "========== Cavity section" << endl;
        os << "Cavity type: GePol" << endl;
	os << "Number of spheres: " << nSpheres << endl;
        os << "Number of finite elements: " << nElements << endl;
/* RDR To be revised...
	for(int i = 0; i < nSpheres + addedSpheres; i++) {
		if ( i < nSpheres ) {
			os << endl;
			os << "Primary Spheres" << endl;
			os << "Sphere " << i+1 << endl;
			os << spheres[i] << endl;
                        os << "Sphere Colour" << endl;
                        os << spheres[i].getSphereColour() << endl;
		} else {
			os << endl;
			os << "Secondary Spheres" << endl;
			os << "Sphere " << i+1 << endl;
			os << spheres[i] << endl;
		}
	}
*/
	return os;
}

void GePolCavity::setMaxAddedSpheres(bool add, int maxAdd) {
	if ( add ) 
	{
		maxAddedSpheres = maxAdd;
	}
	else
	{
		maxAddedSpheres = 0;
	}
}

void GePolCavity::setProbeRadius( double rsolv ) {
	probeRadius = rsolv;
	addSpheres = true;
	maxAddedSpheres = 100;
}
