/*

Wavelet Cavity c++ interface and wrapper methods

*/
#include <iostream>
#include <fstream> 
#include <string>

#include <Eigen/Dense>

using namespace std;

extern "C"{
	//#include "WEM.h"
	//#include "read_points.h"
#include "vector2.h"
#include "vector3.h"
#include "interpolate.h"
#include "interpolate_pwl.h"
//#include "topology.h"
//#include "kern.h"
//#include "compression.h"
//#include "postproc.h"
//#include "WEMRHS.h"
//#include "WEMPCG.h"
//#include "WEMPGMRES.h"
//#include "dwt.h"
#include "cubature.h"
#include "gauss_square.h"
#include "constants.h"
}

#include "PhysicalConstants.hpp"
//#include "Getkw.h"
#include "WaveletCavity.hpp"

/*
WaveletCavity::WaveletCavity(const Getkw & Input, const string path){
	Section cavity = Input.getSect(path);
	vector<double> spheresInput = cavity.getDblVec("Spheres");
	patchLevel = cavity.getInt("PatchLevel");
	probeRadius = cavity.getDbl("ProbeRadius");
	coarsity = cavity.getDbl("Coarsity");
	nSpheres = spheresInput.size()/4; // the correctness of the size has ben checked at input parsing
	Eigen::Matrix3Xd sphereCenter;
	Eigen::VectorXd sphereRadius;
        sphereCenter.resize(Eigen::NoChange, nSpheres);
        sphereRadius.resize(nSpheres);
	int j = 0;
	for (int i = 0; i < nSpheres; i++) {
		sphereCenter(0,i) = spheresInput[j];
		sphereCenter(1,i) = spheresInput[j+1];
		sphereCenter(2,i) = spheresInput[j+2];
		sphereRadius(i)   = spheresInput[j+3];
		j += 4;
	}
	uploadedDyadic = false;
}

WaveletCavity::WaveletCavity(const Section & cavity){
	vector<double> spheresInput = cavity.getDblVec("Spheres");
	patchLevel = cavity.getInt("PatchLevel");
	probeRadius = cavity.getDbl("ProbeRadius");
	coarsity = cavity.getDbl("Coarsity");
	nSpheres = spheresInput.size()/4; // the correctness of the size has ben checked at input parsing
	Eigen::Matrix3Xd sphereCenter;
	Eigen::VectorXd sphereRadius;
        sphereCenter.resize(Eigen::NoChange, nSpheres);
        sphereRadius.resize(nSpheres);
	int j = 0;
	for (int i = 0; i < nSpheres; i++) {
		sphereCenter(0,i) = spheresInput[j];
		sphereCenter(1,i) = spheresInput[j+1];
		sphereCenter(2,i) = spheresInput[j+2];
		sphereRadius(i)   = spheresInput[j+3];
		j += 4;
	}
	uploadedDyadic = false;
}*/


void WaveletCavity::writeInput(string &fileName){

    ofstream output;
    output.open(fileName.c_str(), fstream::out);

    output << nSpheres << endl;
	output.setf(ios_base::showpoint);
	output.precision(12);
    for(int i=0; i < nSpheres; i++) {
		output << sphereCenter(0,i) << " ";
		output << sphereCenter(1,i) << " ";
		output << sphereCenter(2,i) << " ";
		output << sphereRadius(i) << endl;
    }
    output.close();
	uploadedDyadic = false;
}

extern "C" {
	int waveletCavityDrv_(double probeRadius, double coarsity, 
						  int patchLevel, const char* infile);
}

void WaveletCavity::makeCavity() {
	int dummy = 0, check = 0;
	string infile = "cavity.inp";
	writeInput(infile);
	check = waveletCavityDrv_(probeRadius, coarsity, patchLevel, 
							  infile.c_str());
	if (check != 0) {
		std::cout << "Problem with the wavelet cavity!" << std::endl;
		exit(-1);
	}
}


void WaveletCavity::readCavity(const string & filename) {

		int i, j, k;
		double x, y, z;

		ifstream file;
		file.open(filename.c_str());
		file >> nLevels >> nPatches;

		int nNodes = (1 << nLevels) + 1;

		nPoints = nPatches * nNodes * nNodes;
		
		for (int k = 0; k < nPoints; k++) {
			file >> i >> j >> k >> x >> y >> z;
			Eigen::Vector3i index(i, j, k);
			nodeIndex.push_back(index);
			Eigen::Vector3d point(x, y, z);
			nodePoint.push_back(point);
		}

		file.close();
		uploadedDyadic = true;

}

void WaveletCavity::uploadPoints(int quadLevel, vector3 **** T_, bool isPWL) {
	if(isPWL) {
		uploadPointsPWL(quadLevel, T_);
	} else {
		uploadPointsPWC(quadLevel, T_);
	}
}

void WaveletCavity::uploadPointsPWC(int quadLevel, vector3 **** T_) {
	if (not uploadedDyadic) {
		cout << "Error: upload dyadic file first" << endl;
		exit(-1);
	}
	vector2 s, t;
	vector3 point;
	vector3 norm;
	int n = 1 << nLevels;
	double h = 1.0 / n;
	cubature *Q;
	init_Gauss_Square(&Q, quadLevel + 1);

	nElements = nPatches * n * n * Q[quadLevel].nop;

	elementCenter.resize(Eigen::NoChange, nElements);
	elementNormal.resize(Eigen::NoChange, nElements);
	elementArea.resize(nElements);

	int j = 0;
	for (int i1 = 0; i1 < nPatches; i1++){
		for (int i2 = 0; i2 < n; i2++){
			s.y = h * i2;
			for (int i3=0; i3 < n; i3++){
				s.x = h * i3;
				for (int k = 0; k < Q[quadLevel].nop; k++){
					t = vector2_add(s,vector2_Smul(h,Q[quadLevel].xi[k]));
					point = Chi(t,T_[i1], nLevels);
					norm = n_Chi(t,T_[i1], nLevels);
					Eigen::Vector3d center(point.x, point.y, point.z);	 
					Eigen::Vector3d normal(norm.x,  norm.y,  norm.z);	 
					normal.normalize();
					double area = h * h * Q[quadLevel].w[k] * vector3_norm(n_Chi(t, T_[i1], nLevels));
					elementCenter.col(j) = center.transpose();
					elementNormal.col(j) = normal.transpose();
					elementArea(j) = area;
					j++;
				}
			}
		}
	}
	free_Gauss_Square(&Q,quadLevel+1);  
	built = true;
}

void WaveletCavity::uploadPointsPWL(int quadLevel, vector3 **** T_) {
	if (not uploadedDyadic) {
		cout << "Error: upload dyadic file first" << endl;
		exit(-1);
	}
	vector2 s, t;
	vector3 point;
	vector3 norm;
	int n = 1 << nLevels;
	double h = 1.0 / n;
	cubature *Q;
	init_Gauss_Square(&Q, quadLevel + 1);

	nElements = nPatches * n * n * Q[quadLevel].nop;

	elementCenter.resize(Eigen::NoChange, nElements);
	elementNormal.resize(Eigen::NoChange, nElements);
	elementArea.resize(nElements);

	int j = 0;
	for (int i1 = 0; i1 < nPatches; i1++){
		for (int i2 = 0; i2 < n; i2++){
			s.y = h * i2;
			for (int i3=0; i3 < n; i3++){
				s.x = h * i3;
				for (int k = 0; k < Q[quadLevel].nop; k++){
					t = vector2_add(s,vector2_Smul(h,Q[quadLevel].xi[k]));
					point = Chi_pwl(t,T_[i1], nLevels);
					norm = n_Chi_pwl(t,T_[i1], nLevels);
					Eigen::Vector3d center(point.x, point.y, point.z);	 
					Eigen::Vector3d normal(norm.x,  norm.y,  norm.z);	 
					normal.normalize();
					double area = h * h * Q[quadLevel].w[k] * vector3_norm(n_Chi_pwl(t, T_[i1], nLevels));
					elementCenter.col(j) = center.transpose();
					elementNormal.col(j) = normal.transpose();
					elementArea(j) = area;
					j++;
				}
			}
		}
	}
	free_Gauss_Square(&Q,quadLevel+1);  
	built = true;
}

std::ostream & WaveletCavity::printCavity(std::ostream & os) 
{
	os << "========== Cavity section" << endl;
        os << "Cavity type: Wavelet" << endl;
	os << "Probe Radius:   " << probeRadius << endl;
	os << "Coarsity:       " << coarsity << endl;
	os << "Patch Level:    " << patchLevel << endl;
	os << "Number of spheres: " << nSpheres;
        os << "Number of finite elements: " << nElements << endl;
        /*for(int i = 0; i < nElements; i++) 
	{
		os << std::endl;
		os << i+1 << " ";
		os << elementCenter(0,i) << " ";
		os << elementCenter(1,i) << " ";
		os << elementCenter(2,i) << " ";
		os << elementArea(i) << " ";
        }
    	for(int i = 0; i < nSpheres; i++) 
	{
		os << endl;
		os << i+1 << " ";
		os << sphereCenter(0,i) << " ";
		os << sphereCenter(1,i) << " ";
		os << sphereCenter(2,i) << " ";
		os << sphereRadius(i) << " ";
    	}*/
	if (uploadedDyadic) 
	{
		for(int i = 0; i < nPoints; i++) 
		{
			os << endl;
			os << i+1 << " ";
			os << nodeIndex[i].transpose() << " " << nodePoint[i].transpose() << " ";
		}
	}
	return os;
}

