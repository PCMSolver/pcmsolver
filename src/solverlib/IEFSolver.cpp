/*! \file IEFSolver.cpp 
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
#include "IEFSolver.h"

IEFSolver::IEFSolver(GreensFunction & gfi, GreensFunction & gfo) : 
	PCMSolver(gfi, gfo) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver::IEFSolver(GreensFunction * gfi, GreensFunction * gfo) :
	PCMSolver(gfi, gfo) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver::IEFSolver(Section solver) : PCMSolver(solver) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver::~IEFSolver(){
}

void IEFSolver::buildSystemMatrix(Cavity & cavity) {
    if (GePolCavity *gePolCavity = dynamic_cast<GePolCavity*> (&cavity)) {
		if (greenInside->isUniform() && greenOutside->isUniform()) {
			buildIsotropicMatrix(*gePolCavity);
		} else {
			buildAnisotropicMatrix(*gePolCavity);
		}
	} else {
		exit(-1);
	}
}

double IEFSolver::compDiagonalElementSoper(GreensFunction *green, int i, GePolCavity cav) {
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

double IEFSolver::compDiagonalElementDoper(GreensFunction *green, int i, GePolCavity cav) {
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

void IEFSolver::buildAnisotropicMatrix(GePolCavity cav){

    cavitySize = cav.size();

    MatrixXd SI(cavitySize, cavitySize);
    MatrixXd SE(cavitySize, cavitySize);
    MatrixXd DI(cavitySize, cavitySize);
    MatrixXd DE(cavitySize, cavitySize);
    
    for(int i = 0; i < cavitySize; i++){
		Vector3d p1 = cav.getTessCenter(i);
		Vector3d n1 = cav.getTessNormal(i);
		SI(i,i) =  compDiagonalElementSoper(greenInside,  i, cav); 
		SE(i,i) =  compDiagonalElementSoper(greenOutside, i, cav);
		DI(i,i) =  compDiagonalElementDoper(greenInside,  i, cav); 
		DE(i,i) =  compDiagonalElementDoper(greenOutside, i, cav); 
		for (int j = 0; j < cavitySize; j++){
			Vector3d p2 = cav.getTessCenter(j);
			Vector3d n2 = cav.getTessNormal(j);
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
	builtIsotropicMatrix = false;

}

void IEFSolver::buildIsotropicMatrix(GePolCavity cav){


	double epsilon;
    if (UniformDielectric *uniform = 
		dynamic_cast<UniformDielectric *>(greenOutside)) {
		epsilon = uniform->getEpsilon();
	} else {
		cout << "Need uniform dielectric outside" << endl;
		exit(1);
	}
    if (Vacuum *vacuum = dynamic_cast<Vacuum *>(greenInside)) {
	} else {
		cout << "Need vacuum inside" << endl;
		exit(1);
	}

    cavitySize = cav.size();

    MatrixXd SI(cavitySize, cavitySize);
    MatrixXd DI(cavitySize, cavitySize);
    
    for(int i = 0; i < cavitySize; i++){
		Vector3d p1 = cav.getTessCenter(i);
		Vector3d n1 = cav.getTessNormal(i);
		SI(i,i) =  compDiagonalElementSoper(greenInside,  i, cav); 
		DI(i,i) =  compDiagonalElementDoper(greenInside,  i, cav);
		for (int j = 0; j < cavitySize; j++){
			Vector3d p2 = cav.getTessCenter(j);
			Vector3d n2 = cav.getTessNormal(j);
			if (i != j) {
				SI(i,j) = greenInside->evalf(p1, p2);
				DI(i,j) = -greenInside->evald(n2, p1, p2);
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

	double fact = (epsilon+1.0)/(epsilon-1.0);

    PCMMatrix = (fact * aInv - DI) * a * SI;
    PCMMatrix = PCMMatrix.inverse();
    PCMMatrix *= (aInv - DI);
    PCMMatrix = PCMMatrix * a;
	
	builtIsotropicMatrix = true;
	builtAnisotropicMatrix = false;

}

VectorXd IEFSolver::compCharge(const VectorXd &potential) {
	VectorXd charge;
	if (builtIsotropicMatrix or builtAnisotropicMatrix) {
		charge = - PCMMatrix * potential;
	} else {
		cout << "PCM matrix not initialized" << endl;
		exit(1);
	}
	return charge;
}

void IEFSolver::compCharge(const VectorXd & potential, VectorXd & charge) {
	if (builtIsotropicMatrix or builtAnisotropicMatrix) {
		charge = - PCMMatrix * potential;
	} else {
		cout << "PCM matrix not initialized" << endl;
		exit(1);
	}
}

ostream & IEFSolver::printSolver(ostream & os) {
	string type = "Traditional";
	os << "~~~~~~~~~~ PCMSolver ~~~~~~~~~~\n" << endl;
	os << "========== Solver section" << endl;
	os << "Solver Type: " << type << endl;
	os << "Solvent: " << solvent << endl;
	return os;
}
