#include "IEFSolver.hpp"

#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

#include "Cavity.hpp"
#include "GreensFunction.hpp"

void IEFSolver::buildSystemMatrix(Cavity & cavity) 
{
	if (greenInside->isUniform() && greenOutside->isUniform())
	{
		buildIsotropicMatrix(cavity);
	}
	else
	{
		buildAnisotropicMatrix(cavity);
	}
}

void IEFSolver::buildAnisotropicMatrix(Cavity & cav)
{
	int cavitySize = cav.size();
    	Eigen::MatrixXd SI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
	Eigen::MatrixXd SE = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd DI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd DE = Eigen::MatrixXd::Zero(cavitySize, cavitySize);

    	// This is the very core of PCMSolver
    	greenInside->compOffDiagonal(cav.getElementCenter(), cav.getElementNormal(), SI, DI);
    	greenInside->compDiagonal(cav.getElementArea(), cav.getElementRadius(), SI, DI);
    	greenOutside->compOffDiagonal(cav.getElementCenter(), cav.getElementNormal(), SE, DE);
    	greenOutside->compDiagonal(cav.getElementArea(), cav.getElementRadius(), SE, DE);
/*    for(int i = 0; i < cavitySize; i++)
    {
	    Eigen::Vector3d p1 = cav.getElementCenter(i); 
		double area = cav.getElementArea(i);
		double radius = cav.getElementRadius(i);
		SI(i,i) =  greenInside->compDiagonalElementS(area); 
		SE(i,i) = greenOutside->compDiagonalElementS(area); 
		DI(i,i) =  greenInside->compDiagonalElementD(area, radius); 
		DE(i,i) = greenOutside->compDiagonalElementD(area, radius); 
		for (int j = 0; j < cavitySize; j++){
			Eigen::Vector3d p2 = cav.getElementCenter(j);
			Eigen::Vector3d n2 = cav.getElementNormal(j);
			n2.normalize();
			if (i != j) {
				SI(i,j) = greenInside->evalf(p1, p2);
				SE(i,j) = greenOutside->evalf(p1, p2);
				DI(i,j) = greenInside->evald(n2, p1, p2);
				DE(i,j) = greenOutside->evald(n2, p1, p2);
			}
		}
    }*/
    	Eigen::MatrixXd a = cav.getElementArea().asDiagonal();// Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd aInv = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	aInv = 2 * M_PI * a.inverse();
    //a.setZero();
    //aInv.setZero();
/*    for (int i = 0; i < cavitySize; i++) {
		a(i,i) = cav.getElementArea(i);
		aInv(i,i) = 2 * M_PI / cav.getElementArea(i);
    }*/
    	PCMMatrix = ((aInv - DE) * a * SI + SE * a * (aInv + DI.transpose().eval()));
    	PCMMatrix = PCMMatrix.inverse();
    	PCMMatrix *= ((aInv - DE) - SE * SI.inverse() * (aInv - DI));
    	PCMMatrix = PCMMatrix * a;
    	builtAnisotropicMatrix = true;
    	builtIsotropicMatrix = false;
}


void IEFSolver::buildIsotropicMatrix(Cavity & cav)
{
	double epsilon = greenOutside->getDielectricConstant();
        int cavitySize = cav.size();
    	Eigen::MatrixXd SI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd DI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);

    	// This is the very core of PCMSolver
    	greenInside->compOffDiagonal(cav.getElementCenter(), cav.getElementNormal(), SI, DI);
    	greenInside->compDiagonal(cav.getElementArea(), cav.getElementRadius(), SI, DI);

    	Eigen::MatrixXd a = cav.getElementArea().asDiagonal();
    	Eigen::MatrixXd aInv = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	aInv = 2 * M_PI * a.inverse();
	
	double fact = (epsilon + 1.0)/(epsilon - 1.0);
    	PCMMatrix = (fact * aInv - DI) * a * SI;
    	PCMMatrix = PCMMatrix.inverse();
    	PCMMatrix *= (aInv - DI);
    	PCMMatrix = PCMMatrix * a;
    	Eigen::MatrixXd PCMAdjoint(cavitySize, cavitySize); 
    	PCMAdjoint = PCMMatrix.adjoint().eval(); // See Eigen doc for the reason of this
    	PCMMatrix = 0.5 * (PCMMatrix + PCMAdjoint);
	// PRINT TO FILE RELEVANT INFO ABOUT PCMMatrix
    	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(PCMMatrix);
    	if (solver.info() != Eigen::Success) 
		abort();
    	std::ofstream matrixOut("PCM_matrix");
    	matrixOut << "PCM matrix printout" << std::endl;
 	matrixOut << "Number of Tesserae: " << cavitySize << std::endl;
    	matrixOut << "Largest Eigenvalue: " << solver.eigenvalues()[cavitySize-1] << std::endl;
    	matrixOut << "Lowest Eigenvalue: " << solver.eigenvalues()[0] << std::endl;
    	matrixOut << "Average of Eigenvalues: " << (solver.eigenvalues().sum() / cavitySize)<< std::endl;
    	matrixOut << "List of Eigenvalues:\n" << solver.eigenvalues() << std::endl;
    	matrixOut.close();
	builtIsotropicMatrix = true;
	builtAnisotropicMatrix = false;
}

void IEFSolver::compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge) 
{
	if (builtIsotropicMatrix or builtAnisotropicMatrix) 
	{
		charge = - PCMMatrix * potential;
	} 
	else 
	{
		throw std::runtime_error("PCM matrix not initialized!");
	}
}

std::ostream & IEFSolver::printSolver(std::ostream & os) 
{
	std::string type;
	if (builtAnisotropicMatrix) 
	{
		type = "IEFPCM, anisotropic";
	} 
	else 
	{
		type = "IEFPCM, isotropic";
	}
	os << "Solver Type: " << type;
	return os;
}
