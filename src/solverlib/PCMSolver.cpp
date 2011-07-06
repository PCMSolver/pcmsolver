/*! \file PCMSolver.cpp 
\brief PCM solver
*/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "GreensFunctionSum.h"
#include "MetalSphere.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "PCMSolver.h"

GreensFunction& PCMSolver::getGreenInside(){
	return *greenInside;
}

GreensFunction& PCMSolver::getGreenOutside(){
	return *greenOutside;
}

double PCMSolver::compDiagonalElementSoper(GreensFunction *green, int i, GePolCavity cav) {
    double s;
    if (UniformDielectric *uniform = dynamic_cast<UniformDielectric*> (green)) {
	    double eps = uniform->getEpsilon();
	    s = factor * sqrt(4 * M_PI / cav.getTessArea(i)) / eps;   
    }
    else if (Vacuum *vacuum = dynamic_cast<Vacuum *>(green)) {
		s = factor * sqrt(4 * M_PI / cav.getTessArea(i));   
    }
    else {
	    cout << "Not uniform dielectric" << endl;
	    cout << "Not yet implemented" << endl;
	    exit(-1);
    }
    return s;
}

double PCMSolver::compDiagonalElementDoper(GreensFunction *green, int i, GePolCavity cav) {
    double s, d;
    if (UniformDielectric *uniform = dynamic_cast<UniformDielectric *>(green)) {
        s = factor * sqrt(4 * M_PI / cav.getTessArea(i));   
		d = -s / (2*cav.getTessRadius(i));
    }
    else if (Vacuum *vacuum = dynamic_cast<Vacuum *>(green)) {
		s = factor * sqrt(4 * M_PI / cav.getTessArea(i));   
		d = -s / (2*cav.getTessRadius(i));
    }
    else {
		cout << "Not uniform dielectric" << endl;
		cout << "Not yet implemented" << endl;
		exit(-1);
    }
    return d;
}

void PCMSolver::buildPCMMatrix(GePolCavity cav){

    int cavitySize = cav.getNTess();

    MatrixXd SI(cavitySize, cavitySize);
    MatrixXd SE(cavitySize, cavitySize);
    MatrixXd DI(cavitySize, cavitySize);
    MatrixXd DE(cavitySize, cavitySize);
    
    for(int i = 0; i < cavitySize; i++){
		Vector3d p1 = cav.getTessCenter().row(i);
		Vector3d n1 = cav.getTessNormal().row(i);
		SI(i,i) =  compDiagonalElementSoper(greenInside,  i, cav); 
		SE(i,i) =  compDiagonalElementSoper(greenOutside, i, cav);
		DI(i,i) =  compDiagonalElementDoper(greenInside,  i, cav); 
		DE(i,i) =  compDiagonalElementDoper(greenOutside, i, cav); 
		for (int j = 0; j < cavitySize; j++){
			Vector3d p2 = cav.getTessCenter().row(j);
			Vector3d n2 = cav.getTessNormal().row(i);
			if (i != j) {
				SI(i,j) = greenInside->evalf(p1, p2);
				SE(i,j) = greenOutside->evalf(p1, p2);
				DI(i,j) = -greenInside->evald(n2, p1, p2);
				DE(i,j) = -greenOutside->evald(n2, p1, p2);
			}
		}
    }
  
    MatrixXd a(cavitySize, cavitySize);
    MatrixXd aInv(cavitySize, cavitySize);
    a.setZero();
    aInv.setZero();
	cout << cav.getTessArea(1) << endl;

    for (int i = 0; i < cavitySize; i++) {
		a(i,i) = cav.getTessArea(i);
		aInv(i,i) = 2 * M_PI / cav.getTessArea(i);
    }

    PCMMatrix = ((aInv - DE) * a * SI + SE * a * (aInv + DI.transpose()));
    PCMMatrix = PCMMatrix.inverse();
    PCMMatrix *= ((aInv - DE) - SE * SI.inverse() * (aInv - DI));
    PCMMatrix = PCMMatrix * a;
}

bool PCMSolver::readCavity(string &filename){
    Vector3d v;
    Vector3d sph;
    double area;
    double radius;
    int TesNum;
    ifstream input;
    std::cout << "filename " << filename << std::endl;;

    input.open(filename.c_str(), fstream::in);
    if (input.eof()){
	cout << "Unexpected end of file." ;
	exit(1);
    }
    if (input.bad() != 0 ){
	cout << "Type mismatch or file corrupted." ;
	exit(1);
    }


    input >> cavitySize;

    areaTess.resize(cavitySize);
    radiusTess.resize(cavitySize);
    centerSphereTess.resize(cavitySize, 3);
    centerTess.resize(cavitySize, 3);
    PCMMatrix.resize(cavitySize, cavitySize);

    for(int i = 0; i < cavitySize; i++) {
	if (input.eof()){
	    cout << "Unexpected end of file." ;
	    exit(1);
	}
	if (input.bad() != 0 ){
	    cout << "Type mismatch or file corrupted." ;
	    exit(1);
	}
	input >> centerTess(i,0);
	input >> centerTess(i,1);
	input >> centerTess(i,2);
	input >> areaTess(i);
	input >> centerSphereTess(i,0);
	input >> centerSphereTess(i,1);
	input >> centerSphereTess(i,2);
	input >> radiusTess(i);
    }
    input.close();

    normalTess = centerTess - centerSphereTess;
    for(int i = 0; i < cavitySize; i++) {
	normalTess.row(i).normalize();
    }

    std::cout<< centerTess << std::endl;

    return false;
}

//template class PCMSolver<Vacuum, UniformDielectric>;
//template class PCMSolver<Vacuum, GreensFunctionSum>;

