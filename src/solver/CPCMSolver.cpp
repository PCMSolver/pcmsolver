/*! \file CPCMSolver.cpp 
\brief Solver for conductor-like approximation: C-PCM (COSMO)
*/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

#include "Getkw.h"
#include "taylor.hpp"
#include "GreensFunctionInterface.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "PCMSolver.h"
#include "CPCMSolver.h"

CPCMSolver::CPCMSolver(GreensFunctionInterface & gfi, GreensFunctionInterface & gfo) : 
	PCMSolver(gfi, gfo) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

CPCMSolver::CPCMSolver(GreensFunctionInterface * gfi, GreensFunctionInterface * gfo) :
	PCMSolver(gfi, gfo) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

CPCMSolver::CPCMSolver(const Section & solver) : PCMSolver(solver) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

CPCMSolver::~CPCMSolver(){
}

void CPCMSolver::buildSystemMatrix(Cavity & cavity) {
    if (GePolCavity *gePolCavity = dynamic_cast<GePolCavity*> (&cavity)) {
		if (this->greenInside->isUniform() && this->greenOutside->isUniform()) {
			buildIsotropicMatrix(*gePolCavity);
		} else {
			std::cout << "C-PCM is defined only for isotropic environments!" << std::endl;
			exit(-1);
		}
	} else {
		std::cout << "No other cavity than GePol for traditional PCM" << std::endl;
		exit(-1);
	}
}

void CPCMSolver::buildAnisotropicMatrix(GePolCavity & cav){
    cavitySize = cav.size();
    MatrixXd SI(cavitySize, cavitySize);
    MatrixXd SE(cavitySize, cavitySize);
    MatrixXd DI(cavitySize, cavitySize);
    MatrixXd DE(cavitySize, cavitySize);
    for(int i = 0; i < cavitySize; i++){
		Vector3d p1 = cav.getTessCenter(i);
		double area = cav.getTessArea(i);
		double radius = cav.getTessRadius(i);
		SI(i,i) =  greenInside->compDiagonalElementS(area); 
		SE(i,i) = greenOutside->compDiagonalElementS(area); 
		DI(i,i) =  greenInside->compDiagonalElementD(area, radius); 
		DE(i,i) = greenOutside->compDiagonalElementD(area, radius); 
		for (int j = 0; j < cavitySize; j++){
			Vector3d p2 = cav.getTessCenter(j);
			Vector3d n2 = cav.getTessNormal(j);
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


void CPCMSolver::buildIsotropicMatrix(GePolCavity & cav){
    double epsilon = this->greenOutside->getDielectricConstant();
    cavitySize = cav.size();
    Eigen::MatrixXd SI(cavitySize, cavitySize);
    for(int i = 0; i < cavitySize; i++){
		Eigen::Vector3d p1 = cav.getTessCenter(i);
		double area = cav.getTessArea(i);
		double radius = cav.getTessRadius(i);
		SI(i,i) = greenInside->compDiagonalElementS(area); 
		for (int j = 0; j < cavitySize; j++){
			Eigen::Vector3d p2 = cav.getTessCenter(j);
			Eigen::Vector3d n2 = cav.getTessNormal(j);
			n2.normalize();
			if (i != j) {
				SI(i,j) = greenInside->evalf(p1, p2);
			}
		}
    }
    double fact = (epsilon - 1.0)/(epsilon + correction);
    PCMMatrix = SI;
    PCMMatrix = fact * PCMMatrix.inverse();
    Eigen::MatrixXd PCMAdjoint(cavitySize, cavitySize); 
    PCMAdjoint = PCMMatrix.adjoint().eval(); // See Eigen doc for the reason of this
    PCMMatrix = 0.5 * (PCMMatrix + PCMAdjoint);
// PRINT TO FILE RELEVANT INFO ABOUT PCMMatrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(PCMMatrix);
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

void CPCMSolver::compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge) {
	if (builtIsotropicMatrix or builtAnisotropicMatrix) {
		charge = - PCMMatrix * potential;
	} else {
		cout << "PCM matrix not initialized" << endl;
		exit(1);
	}
}
    
void CPCMSolver::setCorrection(double _correction){
	correction = _correction;
}

std::ostream & CPCMSolver::printObject(ostream & os) {
	std::string type = "C-PCM";
	os << "~~~~~~~~~~ PCMSolver ~~~~~~~~~~\n" << endl;
	os << "========== Solver section" << endl;
	os << "Solver Type: " << type << endl;
	os << solvent << endl;
	return os;
}

