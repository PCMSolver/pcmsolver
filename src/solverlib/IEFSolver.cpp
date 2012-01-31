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
#include "GreensFunction.h"
#include "Vacuum.h"
#include "UniformDielectric.h"
#include "GreensFunctionSum.h"
#include "MetalSphere.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "PCMSolver.h"
#include "IEFSolver.h"

IEFSolver<T>::IEFSolver(GreensFunction<T> & gfi, GreensFunction<T> & gfo) : 
	PCMSolver<T>(gfi, gfo) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver<T>::IEFSolver(GreensFunction<T> * gfi, GreensFunction<T> * gfo) :
	PCMSolver<T>(gfi, gfo) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver<T>::IEFSolver(Section solver) : PCMSolver<T>(solver) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver<T>::~IEFSolver(){
}

void IEFSolver<T>::buildSystemMatrix(Cavity & cavity) {
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

double IEFSolver<T>::compDiagonalElementSoper(GreensFunction<T> *green, int i, GePolCavity cav) {
    double s;
    if (UniformDielectric<T> * uniform = dynamic_cast<UniformDielectric<T> * > (green)) {
	    double eps = uniform->getEpsilon();
		s = factor * sqrt(4 * M_PI / (cav.getTessArea)(i)) / eps;   
    } else if (Vacuum<T> *vacuum = dynamic_cast<Vacuum<T> * >(green) ) {
		s = factor * sqrt(4 * M_PI / (cav.getTessArea)(i));   
    } else {
	    cout << "Not uniform dielectric" << endl;
	    cout << "Not yet implemented" << endl;
	    exit(-1);
    }
    return s;
}

double IEFSolver<T>::compDiagonalElementDoper(GreensFunction<T> *green, int i, GePolCavity cav) {
    double s, d;
    if (UniformDielectric<T> *uniform = dynamic_cast<UniformDielectric<T> *>(green)) {
        s = factor * sqrt(4 * M_PI / cav.getTessArea(i));   
		d = -s / (2*cav.getTessRadius(i));
    } else if (Vacuum<T>* vacuum = dynamic_cast<Vacuum<T> *>(green)) {
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

void IEFSolver<T>::buildAnisotropicMatrix(GePolCavity cav){

    this->cavitySize = cav.size();

    MatrixXd SI(this->cavitySize, this->cavitySize);
    MatrixXd SE(this->cavitySize, this->cavitySize);
    MatrixXd DI(this->cavitySize, this->cavitySize);
    MatrixXd DE(this->cavitySize, this->cavitySize);
    
    for(int i = 0; i < this->cavitySize; i++){
		Vector3d p1 = cav.getTessCenter(i);
		Vector3d n1 = cav.getTessNormal(i);
		SI(i,i) =  compDiagonalElementSoper(this->greenInside,  i, cav); 
		SE(i,i) =  compDiagonalElementSoper(this->greenOutside, i, cav);
		DI(i,i) =  compDiagonalElementDoper(this->greenInside,  i, cav); 
		DE(i,i) =  compDiagonalElementDoper(this->greenOutside, i, cav); 
		for (int j = 0; j < this->cavitySize; j++){
			Vector3d p2 = cav.getTessCenter(j);
			Vector3d n2 = cav.getTessNormal(j);
			if (i != j) {
				SI(i,j) =   this->greenInside->evalf(p1, p2);
				SE(i,j) =   this->greenOutside->evalf(p1, p2);
				DI(i,j) = - this->greenInside->evald(n2, p1, p2);
				DE(i,j) = - this->greenOutside->evald(n2, p1, p2);
			}
		}
    }
  
    MatrixXd a(this->cavitySize, this->cavitySize);
    MatrixXd aInv(this->cavitySize, this->cavitySize);
    a.setZero();
    aInv.setZero();

    for (int i = 0; i < this->cavitySize; i++) {
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

void IEFSolver<T>::buildIsotropicMatrix(GePolCavity cav){
	double epsilon;
    if (UniformDielectric<T> *uniform = 
		dynamic_cast<UniformDielectric<T> *>(this->greenOutside)) {
		epsilon = uniform->getEpsilon();
	} else {
		cout << "Need uniform dielectric outside" << endl;
		exit(1);
	}
	if (Vacuum<T> *vacuum = dynamic_cast<Vacuum<T> *>(this->greenInside)) {
	} else {
		cout << "Need vacuum inside" << endl;
		exit(1);
	}

    this->cavitySize = cav.size();

    MatrixXd SI(this->cavitySize, this->cavitySize);
    MatrixXd DI(this->cavitySize, this->cavitySize);
    
    for(int i = 0; i < this->cavitySize; i++){
		Vector3d p1 = cav.getTessCenter(i);
		Vector3d n1 = cav.getTessNormal(i);
		SI(i,i) =  compDiagonalElementSoper(this->greenInside,  i, cav); 
		DI(i,i) =  compDiagonalElementDoper(this->greenInside,  i, cav);
		for (int j = 0; j < this->cavitySize; j++){
			Vector3d p2 = cav.getTessCenter(j);
			Vector3d n2 = cav.getTessNormal(j);
			if (i != j) {
				SI(i,j) = this->greenInside->evalf(p1, p2);
				DI(i,j) = -this->greenInside->evald(n2, p1, p2);
			}
		}
    }
  
    MatrixXd a(this->cavitySize, this->cavitySize);
    MatrixXd aInv(this->cavitySize, this->cavitySize);
    a.setZero();
    aInv.setZero();

    for (int i = 0; i < this->cavitySize; i++) {
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

template<class T>
VectorXd IEFSolver<T>::compCharge(const VectorXd &potential) {
	VectorXd charge;
	if (builtIsotropicMatrix or builtAnisotropicMatrix) {
		charge = - PCMMatrix * potential;
	} else {
		cout << "PCM matrix not initialized" << endl;
		exit(1);
	}
	return charge;
}

template<class T>
void IEFSolver<T>::compCharge(const VectorXd & potential, VectorXd & charge) {
	if (builtIsotropicMatrix or builtAnisotropicMatrix) {
		charge = - PCMMatrix * potential;
	} else {
		cout << "PCM matrix not initialized" << endl;
		exit(1);
	}
}

template<class T>
ostream & IEFSolver<T>::printObject(ostream & os) {
	string type = "Traditional";
	os << "~~~~~~~~~~ PCMSolver ~~~~~~~~~~\n" << endl;
	os << "========== Solver section" << endl;
	os << "Solver Type: " << type << endl;
	os << "Solvent: " << this->solvent << endl;
	return os;
}

template class IEFSolver<double>;
template class IEFSolver<taylor<double, 1, 1> >;
template class IEFSolver<taylor<double, 3, 1> >;
template class IEFSolver<taylor<double, 3 ,2> >;
