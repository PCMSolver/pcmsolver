/*

Wavelet Cavity c++ interface and wrapper methods

*/
#include <iostream>
#include <fstream> 
#include <string>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "Cavity.h"
#include "WaveletCavity.h"

WaveletCavity::WaveletCavity(const Getkw & Input, const string path){
	Section cavity = Input.getSect(path);
	vector<double> spheresInput = cavity.getDblVec("Spheres");
	patchLevel = cavity.getInt("PatchLevel");
	probeRadius = cavity.getDbl("ProbeRadius");
	coarsity = cavity.getDbl("Coarsity");
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

WaveletCavity::WaveletCavity(const Section & cavity){
	vector<double> spheresInput = cavity.getDblVec("Spheres");
	patchLevel = cavity.getInt("PatchLevel");
	probeRadius = cavity.getDbl("ProbeRadius");
	coarsity = cavity.getDbl("Coarsity");
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

void WaveletCavity::writeInput(string &fileName){

    ofstream output;
    output.open(fileName.c_str(), fstream::out);

    output << nSpheres << endl;
    for(int i=0; i < nSpheres; i++) {
		output << sphereCenter(0,i) << " ";
		output << sphereCenter(1,i) << " ";
		output << sphereCenter(2,i) << " ";
		output << sphereRadius(i) << endl;
    }
    output.close();
}

extern "C" {
	int waveletCavityDrv_(double probeRadius, double coarsity, 
						  int patchLevel);

}


void WaveletCavity::makeCavity() {
	int dummy = 0, check = 0;
	string fileName = "cavity.inp";
	cout << fileName << endl;
	writeInput(fileName);
	check = waveletCavityDrv_(probeRadius, coarsity, patchLevel);
	cout << "Created wavelet cavity: " << check << endl;
}

ostream & operator<<(ostream &os, const WaveletCavity &cavity) {
	os << "Molecular cavity" << endl;
	os << "Probe Radius:   " << cavity.probeRadius << endl;
	os << "Coarsity:       " << cavity.coarsity << endl;
	os << "Patch Level:    " << cavity.patchLevel << endl;
	os << "Nr. of spheres: " << cavity.nSpheres;
    for(int i = 0; i < cavity.nSpheres; i++) {
		os << endl;
		os << i+1 << " ";
		os << cavity.sphereCenter(0,i) << " ";
		os << cavity.sphereCenter(1,i) << " ";
		os << cavity.sphereCenter(2,i) << " ";
		os << cavity.sphereRadius(i) << " ";
    }
	return os;
}

