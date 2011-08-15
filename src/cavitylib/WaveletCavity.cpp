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

WaveletCavity::WaveletCavity(Getkw &Input){
	vector<double> spheresInput = Input.getDblVec("Cavity.Spheres");
	patchLevel = Input.getInt("Cavity.PatchLevel");
	probeRadius = Input.getDbl("Cavity.ProbeRadius");
	coarsity = Input.getDbl("Cavity.Coarsity");
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
    for(int i=0; i < nTess; i++) {
		output << sphereCenter(0,i) << " ";
		output << sphereCenter(1,i) << " ";
		output << sphereCenter(2,i) << " ";
		output << sphereRadius(i) << endl;
    }
    output.close();
}

void WaveletCavity::makeCavity() {
	int dummy = 0, check = 0;
	string fileName = "cavity.inp";
	writeInput(fileName);
	check = bihp_test_romg(fileName.c_str(), probeRadius, coarsity, &dummy);
	if (check == 1) {
		return 0;
	} else {
		cout << "Error in creating wavelet cavity" << endl; 
		exit(1);
	}
}
