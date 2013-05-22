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
#include "taylor.hpp"
#include "GreensFunctionInterface.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "PCMSolver.h"
#include "IEFSolver.h"

IEFSolver::IEFSolver(GreensFunctionInterface & gfi, GreensFunctionInterface & gfo) : 
	PCMSolver(gfi, gfo) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver::IEFSolver(GreensFunctionInterface * gfi, GreensFunctionInterface * gfo) :
	PCMSolver(gfi, gfo) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver::IEFSolver(const Section & solver) : PCMSolver(solver) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver::~IEFSolver(){
}

void IEFSolver::buildSystemMatrix(Cavity & cavity) {
    if (GePolCavity *gePolCavity = dynamic_cast<GePolCavity*> (&cavity)) {
		if (this->greenInside->isUniform() && this->greenOutside->isUniform()) {
			buildIsotropicMatrix(*gePolCavity);
		} else {
			buildAnisotropicMatrix(*gePolCavity);
		}
	} else {
		std::cout << "No other cavity than GePol for traditional PCM" << std::endl;
		exit(-1);
	}
}

void IEFSolver::buildAnisotropicMatrix(GePolCavity & cav){
    cavitySize = cav.size();
    MatrixXd SI(cavitySize, cavitySize);
    MatrixXd SE(cavitySize, cavitySize);
    MatrixXd DI(cavitySize, cavitySize);
    MatrixXd DE(cavitySize, cavitySize);
    for(int i = 0; i < cavitySize; i++){
		Vector3d p1 = cav.getElementCenter(i);
		double area = cav.getElementArea(i);
		double radius = cav.getElementRadius(i);
		SI(i,i) =  greenInside->compDiagonalElementS(area); 
		SE(i,i) = greenOutside->compDiagonalElementS(area); 
		DI(i,i) =  greenInside->compDiagonalElementD(area, radius); 
		DE(i,i) = greenOutside->compDiagonalElementD(area, radius); 
		for (int j = 0; j < cavitySize; j++){
			Vector3d p2 = cav.getElementCenter(j);
			Vector3d n2 = cav.getElementNormal(j);
			n2.normalize();
			if (i != j) {
				SI(i,j) = greenInside->evalf(p1, p2);
				SE(i,j) = greenOutside->evalf(p1, p2);
				DI(i,j) = greenInside->evald(n2, p1, p2);
				DE(i,j) = greenOutside->evald(n2, p1, p2);
			}
		}
    }
    MatrixXd a(cavitySize, cavitySize);
    MatrixXd aInv(cavitySize, cavitySize);
    a.setZero();
    aInv.setZero();
    for (int i = 0; i < cavitySize; i++) {
		a(i,i) = cav.getElementArea(i);
		aInv(i,i) = 2 * M_PI / cav.getElementArea(i);
    }
    PCMMatrix = ((aInv - DE) * a * SI + SE * a * (aInv + DI.transpose()));
    PCMMatrix = PCMMatrix.inverse();
    PCMMatrix *= ((aInv - DE) - SE * SI.inverse() * (aInv - DI));
    PCMMatrix = PCMMatrix * a;
	builtAnisotropicMatrix = true;
	builtIsotropicMatrix = false;
}


void IEFSolver::buildIsotropicMatrix(GePolCavity & cav){
	double epsilon = this->greenOutside->getDielectricConstant();
    cavitySize = cav.size();
    MatrixXd SI(cavitySize, cavitySize);
    MatrixXd DI(cavitySize, cavitySize);
    for(int i = 0; i < cavitySize; i++){
		Vector3d p1 = cav.getElementCenter(i);
		double area = cav.getElementArea(i);
		double radius = cav.getElementRadius(i);
		SI(i,i) = greenInside->compDiagonalElementS(area); 
		DI(i,i) = greenInside->compDiagonalElementD(area, radius); 
		for (int j = 0; j < cavitySize; j++){
			Vector3d p2 = cav.getElementCenter(j);
			Vector3d n2 = cav.getElementNormal(j);
			n2.normalize();
			if (i != j) {
				SI(i,j) = greenInside->evalf(p1, p2);
				DI(i,j) = greenInside->evald(n2, p1, p2);
			}
		}
    }
    MatrixXd a(cavitySize, cavitySize);
    MatrixXd aInv(cavitySize, cavitySize);
    a.setZero();
    aInv.setZero();
    for (int i = 0; i < cavitySize; i++) {
		a(i,i) = cav.getElementArea(i);
		aInv(i,i) = 2 * M_PI / cav.getElementArea(i);
    }
	double fact = (epsilon+1.0)/(epsilon-1.0);
    PCMMatrix = (fact * aInv - DI) * a * SI;
    PCMMatrix = PCMMatrix.inverse();
    PCMMatrix *= (aInv - DI);
    PCMMatrix = PCMMatrix * a;
    MatrixXd PCMAdjoint(cavitySize, cavitySize); 
    PCMAdjoint = PCMMatrix.adjoint().eval(); // See Eigen doc for the reason of this
    PCMMatrix = 0.5 * (PCMMatrix + PCMAdjoint);
// PRINT TO FILE RELEVANT INFO ABOUT PCMMatrix
    SelfAdjointEigenSolver<MatrixXd> solver(PCMMatrix);
    if (solver.info() != Success) abort();
    ofstream matrixOut("PCM_matrix");
    matrixOut << "PCM matrix printout" << endl;
    matrixOut << "Number of Tesserae: " << cavitySize << endl;
    matrixOut << "Largest Eigenvalue: " << solver.eigenvalues()[cavitySize-1] << endl;
    matrixOut << "Lowest Eigenvalue: " << solver.eigenvalues()[0] << endl;
    matrixOut << "Average of Eigenvalues: " << (solver.eigenvalues().sum() / cavitySize)<< endl;
    matrixOut << "List of Eigenvalues:\n" << solver.eigenvalues() << endl;
    matrixOut.close();
	builtIsotropicMatrix = true;
	builtAnisotropicMatrix = false;
}

void IEFSolver::compCharge(const VectorXd & potential, VectorXd & charge) {
	if (builtIsotropicMatrix or builtAnisotropicMatrix) {
		charge = - PCMMatrix * potential;
	} else {
		cout << "PCM matrix not initialized" << endl;
		exit(1);
	}
}

ostream & IEFSolver::printObject(ostream & os) {
	std::string type;
	if (builtAnisotropicMatrix) {
		type = "IEFPCM, anisotropic";
	} else {
		type = "IEFPCM, isotropic";
	}
	os << "~~~~~~~~~~ PCMSolver ~~~~~~~~~~\n" << endl;
	os << "========== Solver section" << endl;
	os << "Solver Type: " << type << endl;
	os << solvent << endl;
	return os;
}

