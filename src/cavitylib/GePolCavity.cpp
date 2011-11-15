/*

GePol c++ interface and wrapper methods
written by Krzysztof Mozgawa, 2011

*/
#include <iostream>
#include <fstream> 
#include <string>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "Cavity.h"
#include "GePolCavity.h"

GePolCavity::GePolCavity(Getkw &Input){
	vector<double> spheresInput = Input.getDblVec("Cavity.Spheres");
	averageArea = Input.getDbl("Cavity.Area");
	nSpheres = spheresInput.size()/4; // the correctness of the size has ben checked at input parsing
    sphereCenter.resize(NoChange, nSpheres);
    sphereRadius.resize(nSpheres);
	int j = 0;
	for (int i = 0; i < nSpheres; i++) {
		sphereCenter(0,i) = spheresInput[j];
		sphereCenter(1,i) = spheresInput[j+1];
		sphereCenter(2,i) = spheresInput[j+2];
		sphereRadius(i)   = spheresInput[j+3];
		j += 4;
	}
}

bool GePolCavity::readInput(string &filename){

    ifstream input;
    input.open(filename.c_str(), fstream::in);
    if (input.eof()){
		cout << "Unexpected end of file." ;
		exit(1);
    }
    if (input.bad() != 0 ){
		cout << "Type mismatch or file corrupted." ;
		exit(1);
    }
	
    input >> nSpheres;
    sphereCenter.resize(nSpheres, NoChange);
    sphereRadius.resize(nSpheres);
    for(int i = 0; i < nSpheres; i++) {
		if (input.eof()){
			cout << "Unexpected end of file." ;
			exit(1);
		}
		if (input.bad() != 0 ){
			cout << "Type mismatch or file corrupted." ;
			exit(1);
		}
    }
    input.close();
    return false;
}

void GePolCavity::writeOutput(string &filename){

    ofstream output;
    output.open(filename.c_str(), fstream::out);

    output << nTess << endl;
    for(int i=0; i < nTess; i++) {
		output << tessCenter(0,i) << " ";
		output << tessCenter(1,i) << " ";
		output << tessCenter(2,i) << " ";
		output << tessArea(i) << " ";
		output << tessSphereCenter(0,i) << " ";
		output << tessSphereCenter(1,i) << " ";
		output << tessSphereCenter(2,i) << " ";
		output << tessRadius(i) << endl;
    }
    output.close();
}

extern"C" {
    void generatecavity_cpp_(double *xtscor, double *ytscor, double *ztscor, 
							 double *ar, double *xsphcor, double *ysphcor, 
							 double *zsphcor, double *rsph, int *nts, int *nesfp,
							 double *xe, double *ye, double *ze, double *rin, 
							 double *avgArea, double* work, int* lwork);
}


void GePolCavity::makeCavity(int maxts, int lwork){

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

	double *xe = sphereCenter.row(0).data();
	double *ye = sphereCenter.row(1).data();
	double *ze = sphereCenter.row(2).data();
	double *rin = sphereRadius.data();

	generatecavity_cpp_(xtscor, ytscor, ztscor, ar, xsphcor, ysphcor, zsphcor, rsph, &nts, &nSpheres, 
						xe, ye, ze, rin, &averageArea, work, &lwork);

    nTess = int(nts);
    tessCenter.resize(NoChange, nTess);
    tessSphereCenter.resize(NoChange, nTess);
    tessNormal.resize(NoChange, nTess);
    tessArea.resize(nTess);
    tessRadius.resize(nTess);
    for(int i=0; i < nTess; i++){
		tessCenter(0,i) = xtscor[i];
		tessCenter(1,i) = ytscor[i];
		tessCenter(2,i) = ztscor[i];
		tessArea(i) = ar[i];
		tessSphereCenter(0,i) = xsphcor[i];
		tessSphereCenter(1,i) = ysphcor[i];
		tessSphereCenter(2,i) = zsphcor[i];
		tessRadius(i) = rsph[i];
    }

    tessNormal = tessCenter - tessSphereCenter;
    for(int i=0; i < nTess; i++){
		tessNormal.col(i) /= tessNormal.col(i).norm();
	}

	delete xtscor;
	delete ytscor;
	delete ztscor;
	delete ar;
	delete xsphcor;
	delete ysphcor;
	delete zsphcor;
	delete rsph;
	delete work;
	
	isBuilt = true;

}

