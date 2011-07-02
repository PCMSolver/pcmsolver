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
	nSpheres = spheresInput.size()/4; // the correctness of the size has ben checked at input parsing
    sphereCenter.resize(nSpheres, NoChange);
    sphereRadius.resize(nSpheres);
	int j = 0;
	for (int i = 0; i < nSpheres; i++) {
		sphereCenter(i,0) = spheresInput[j];
		sphereCenter(i,1) = spheresInput[j+1];
		sphereCenter(i,2) = spheresInput[j+2];
		sphereRadius(i)   = spheresInput[j+3];
		j += 4;
	}

	//	cout << x << endl;
	//cout << y << endl;
	//cout << z << endl;

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
//		sphereCenter(i,0) = xe[i];
//		sphereCenter(i,1) = ye[i];
//		sphereCenter(i,2) = ze[i];
//		sphereRadius(i) = rin[i];
    }
    input.close();
    return false;
}

void GePolCavity::writeOutput(string &filename){

    ofstream output;
    output.open(filename.c_str(), fstream::out);

    output << nTess << endl;
    for(int i=0; i < nTess; i++) {
		output << tessCenter(i,0) << " ";
		output << tessCenter(i,1) << " ";
		output << tessCenter(i,2) << " ";
		output << tessArea(i) << " ";
		output << tessSphereCenter(i,0) << " ";
		output << tessSphereCenter(i,1) << " ";
		output << tessSphereCenter(i,2) << " ";
		output << tessRadius(i) << endl;
    }
    output.close();
}

extern"C" {
    void generatecavity_cpp_(double *xtscor, double *ytscor, double *ztscor, double *ar, double *xsphcor,
			     double *ysphcor, double *zsphcor, double *rsph, int *nts, int *nesfp,
			     double *xe, double *ye, double *ze, double *rin);
}


void GePolCavity::makeCavity(){

    double xtscor[5000], ytscor[5000], ztscor[5000];
    double xsphcor[5000], ysphcor[5000], zsphcor[5000];
    double ar[5000], rsph[5000];
    int nts;

	double *xe = sphereCenter.col(0).data();
	double *ye = sphereCenter.col(1).data();
	double *ze = sphereCenter.col(2).data();
	double *rin = sphereRadius.data();

	cout << nSpheres << endl;

	for (int i = 0; i < nSpheres; i++) {
		cout << xe[i] << " "  << ye[i] << " "  << ze[i] << " "  << rin[i] << endl;
	}

	generatecavity_cpp_(xtscor, ytscor, ztscor, ar, xsphcor, ysphcor, zsphcor, rsph, &nts, &nSpheres, xe, ye, ze, rin);

	cout << nts << endl;
    
    nTess = int(nts);
    tessCenter.resize(nTess, NoChange);
    tessSphereCenter.resize(nTess,NoChange);
    tessNormal.resize(nTess,NoChange);
    tessArea.resize(nTess);
    tessRadius.resize(nTess);
    for(int i=0; i < nTess; i++){
		tessCenter(i,0) = xtscor[i];
		tessCenter(i,1) = ytscor[i];
		tessCenter(i,2) = ztscor[i];
		tessArea(i) = ar[i];
		tessSphereCenter(i,0) = xsphcor[i];
		tessSphereCenter(i,1) = ysphcor[i];
		tessSphereCenter(i,2) = zsphcor[i];
		tessRadius(i) = rsph[i];
    }
    tessNormal = tessCenter - tessSphereCenter;
}


