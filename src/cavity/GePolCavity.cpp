/*

  GePol c++ interface and wrapper methods
  written by Krzysztof Mozgawa, 2011

*/
#include <iostream>
#include <fstream> 
#include <string>
#include <vector>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "CavityOfSpheres.h"
#include "GePolCavity.h"

GePolCavity::GePolCavity(const Section & cavity){
	averageArea = cavity.getDbl("Area");
	std::string modeString = cavity.getStr("Mode");
	setMode(modeString);
	if (mode == Explicit) {
		vector<double> spheresInput = cavity.getDblVec("Spheres");
		nSpheres = spheresInput.size()/4; // the correctness of the size has ben checked at input parsing
		int j = 0;
		for (int i = 0; i < nSpheres; i++) {
			Vector3d center; 
			center << spheresInput[j], spheresInput[j+1], spheresInput[j+2];
			Sphere sph(center, spheresInput[j+3]);
			spheres.push_back(sph);
			j += 4;
		}
	}
}

GePolCavity::GePolCavity(const Section &cavity, const vector<Sphere> & _spheres){
	averageArea = cavity.getDbl("Area");
	nSpheres = _spheres.size();	
	spheres = _spheres;
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
	
   	VectorXd xv(nSpheres + maxAddedSpheres);
	VectorXd yv(nSpheres + maxAddedSpheres);
	VectorXd zv(nSpheres + maxAddedSpheres);
	VectorXd sphereRadius(nSpheres + maxAddedSpheres);
	
	for ( int i = 0; i < nSpheres; i++ ) {
		xv(i) = spheres[i].getSphereCenter(0);
		yv(i) = spheres[i].getSphereCenter(1);
		zv(i) = spheres[i].getSphereCenter(2);
		sphereRadius(i) = spheres[i].getSphereRadius();
	}
		
	double *xe = xv.data();
	double *ye = yv.data();
	double *ze = zv.data();

	double *rin = sphereRadius.data();
 
	generatecavity_cpp_(xtscor, ytscor, ztscor, ar, xsphcor, ysphcor, zsphcor, rsph, &nts, &nSpheres, 
						xe, ye, ze, rin, &averageArea, &probeRadius, work, &lwork);
	
	VectorXd rtmp(nSpheres + maxAddedSpheres);
	VectorXd xtmp(nSpheres + maxAddedSpheres);
	VectorXd ytmp(nSpheres + maxAddedSpheres);
	VectorXd ztmp(nSpheres + maxAddedSpheres);
	
	rtmp = sphereRadius.head(nSpheres);
	xtmp = xv.head(nSpheres);
	ytmp = yv.head(nSpheres);
	ztmp = zv.head(nSpheres);
	addedSpheres = 0;
	for ( int i = nSpheres; i < nSpheres + addedSpheres; i++ ) {
		if ( sphereRadius(i) != 0) {
			rtmp << sphereRadius(i);
			xtmp << xv(i);
			ytmp << yv(i);
			ztmp << zv(i);
			addedSpheres += 1;
		}
	}
	rtmp.resize(nSpheres + addedSpheres);
	xtmp.resize(nSpheres + addedSpheres);
	ytmp.resize(nSpheres + addedSpheres);
	ztmp.resize(nSpheres + addedSpheres);
	sphereRadius.resize(nSpheres + addedSpheres);
	xv.resize(nSpheres + addedSpheres);
	yv.resize(nSpheres + addedSpheres);
	zv.resize(nSpheres + addedSpheres);
	sphereRadius = rtmp;
    	xv = xtmp;
	yv = ytmp;
	zv = ztmp;
	
	for ( int i = 0; i < nSpheres + addedSpheres; i++ ) {
		Vector3d coord; 
		coord << xv(i), yv(i), zv(i);
		spheres[i].setSphereCenter(coord);
		spheres[i].setSphereRadius(sphereRadius(i));
	}
    	
    nElements = int(nts);
    elementCenter.resize(NoChange, nElements);
    elementSphereCenter.resize(NoChange, nElements);
    elementNormal.resize(NoChange, nElements);
    elementArea.resize(nElements);
    elementRadius.resize(nElements);
    for(int i=0; i < nElements; i++){
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
    for(int i=0; i < nElements; i++){
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
/*
void GePolCavity::setMode(const string & type) {
	if (type == "Atoms") {
		setMode(Atoms);
	} else if (type == "Implicit") {
		setMode(Implicit);
	} else if (type == "Explicit") {
		setMode(Explicit);
	} else {
		exit(-1);
	}
}

void GePolCavity::setMode(int type) {
	switch (type) {
	case Atoms :
		mode = Atoms;
		break;
	case Implicit :
		mode = Implicit;
		break;
	case Explicit :
		mode = Explicit;
		break;
	default :
		exit(-1);
	}
}
*/
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
	if ( add == true ) {
		maxAddedSpheres = maxAdd;
	}
	addSpheres = true;
}

void GePolCavity::setProbeRadius( double rsolv ) {
	probeRadius = rsolv;
	addSpheres = true;
	maxAddedSpheres = 100;
}
