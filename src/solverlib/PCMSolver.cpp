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

#include "Getkw.h"
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "GreensFunctionSum.h"
#include "MetalSphere.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "PCMSolver.h"

PCMSolver::PCMSolver(GreensFunction &gfi, GreensFunction &gfo){
	allocated = false;
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
	greenInside = &gfi; 
	greenOutside = &gfo;
}
PCMSolver::PCMSolver(GreensFunction *gfi, GreensFunction *gfo){
	allocated = false;
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
	greenInside = gfi; 
	greenOutside = gfo;
}

PCMSolver::PCMSolver(Section solver) {
	allocated = true;
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
	greenInside  = 
		greenInside->allocateGreensFunction(solver.getSect("Green<inside>"));
	greenOutside = 
		greenOutside->allocateGreensFunction(solver.getSect("Green<outside>"));
}

PCMSolver::~PCMSolver(){
	if(allocated) {
		delete greenInside; 
		delete greenOutside;
	}
}

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
	    s = factor * sqrt(4 * M_PI / (cav.getTessArea)(i)) / eps;   
    }
    else if (Vacuum *vacuum = dynamic_cast<Vacuum *>(green)) {
		s = factor * sqrt(4 * M_PI / (cav.getTessArea)(i));   
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

void PCMSolver::buildAnisotropicMatrix(GePolCavity cav){

    cavitySize = cav.size();

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

    for (int i = 0; i < cavitySize; i++) {
		a(i,i) = cav.getTessArea(i);
		aInv(i,i) = 2 * M_PI / cav.getTessArea(i);
    }

    PCMMatrix = ((aInv - DE) * a * SI + SE * a * (aInv + DI.transpose()));
    PCMMatrix = PCMMatrix.inverse();
    PCMMatrix *= ((aInv - DE) - SE * SI.inverse() * (aInv - DI));
    PCMMatrix = PCMMatrix * a;
	
	builtAnisotropicMatrix = true;

}

void PCMSolver::buildIsotropicMatrix(GePolCavity cav){
	cout << "Not yet implemented" << endl;
	exit(1);
}

VectorXd PCMSolver::compCharge(const VectorXd &potential) {
	VectorXd charge;
	if (builtIsotropicMatrix or builtAnisotropicMatrix) {
		charge = PCMMatrix * potential;
	} else {
		cout << "PCM matrtix not initialized" << endl;
		exit(1);
	}
	return charge;
}


