#include "IEFSolver.hpp"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
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
#include "MathUtils.hpp"

void IEFSolver::buildSystemMatrix(const Cavity & cavity) 
{
	if (greenInside_->isUniform() && greenOutside_->isUniform())
	{
		buildIsotropicMatrix(cavity);
	}
	else
	{
		buildAnisotropicMatrix(cavity);
	}
}

void IEFSolver::buildAnisotropicMatrix(const Cavity & cav)
{
	int cavitySize = cav.size();
    	Eigen::MatrixXd SI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
	Eigen::MatrixXd SE = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd DI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd DE = Eigen::MatrixXd::Zero(cavitySize, cavitySize);

    	// This is the very core of PCMSolver
    	greenInside_->compOffDiagonal(cav.getElementCenter(), cav.getElementNormal(), SI, DI);
    	greenInside_->compDiagonal(cav.getElementArea(), cav.getElementRadius(), SI, DI);
    	greenOutside_->compOffDiagonal(cav.getElementCenter(), cav.getElementNormal(), SE, DE);
    	greenOutside_->compDiagonal(cav.getElementArea(), cav.getElementRadius(), SE, DE);
    	
	Eigen::MatrixXd a = cav.getElementArea().asDiagonal();// Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd aInv = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	aInv = 2 * M_PI * a.inverse();

    	PCMMatrix = ((aInv - DE) * a * SI + SE * a * (aInv + DI.transpose().eval()));
    	PCMMatrix = PCMMatrix.inverse();
    	PCMMatrix *= ((aInv - DE) - SE * SI.inverse() * (aInv - DI));
    	PCMMatrix = PCMMatrix * a;
    	builtAnisotropicMatrix = true;
    	builtIsotropicMatrix = false;
}


void IEFSolver::buildIsotropicMatrix(const Cavity & cav)
{
	bool symmetry = cav.pointGroup().groupInteger();
	double epsilon = greenOutside_->getDielectricConstant();
        int cavitySize = cav.size();
    	Eigen::MatrixXd SI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	Eigen::MatrixXd DI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);

	// Compute SI and DI on the whole cavity, regardless of symmetry
    	greenInside_->compOffDiagonal(cav.getElementCenter(), cav.getElementNormal(), SI, DI);
    	greenInside_->compDiagonal(cav.getElementArea(), cav.getElementRadius(), SI, DI);
	// Perform symmetry blocking
	// If the group is C1 avoid symmetry blocking, we will just pack the PCMMatrix
	// into "block diagonal" when all other manipulations are done.
	if (symmetry)
	{
		symmetryBlocking(DI, cav);
		symmetryBlocking(SI, cav);
	}

    	Eigen::MatrixXd a = cav.getElementArea().asDiagonal();
    	Eigen::MatrixXd aInv = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	aInv = a.inverse();

	// Tq = -Rv -> q = -(T^-1 * R)v = -Kv
	// T = (2 * M_PI * fact * aInv - DI) * a * SI; R = (2 * M_PI * aInv - DI)
	// PCMMatrix = K = T^-1 * R * a
	// 1. Form T
	double fact = (epsilon + 1.0)/(epsilon - 1.0);
    	PCMMatrix = (2 * M_PI * fact * aInv - DI) * a * SI;
	// 2. Invert T using LU decomposition with full pivoting
	//    This is a rank-revealing LU decomposition, this allows us
	//    to test if T is invertible before attempting to invert it.
	Eigen::FullPivLU<Eigen::MatrixXd> T_LU(PCMMatrix);
	if (!(T_LU.isInvertible()))
		throw std::runtime_error("T matrix is not invertible!");
    	PCMMatrix = T_LU.inverse();
	// 3. Multiply T^-1 and R
    	PCMMatrix *= (2 * M_PI * aInv - DI);
	// 4. Multiply by a
    	PCMMatrix *= a;
	// 5. Symmetrize K := (K + K+)/2
    	Eigen::MatrixXd PCMAdjoint(cavitySize, cavitySize); 
    	PCMAdjoint = PCMMatrix.adjoint().eval(); // See Eigen doc for the reason of this
    	PCMMatrix = 0.5 * (PCMMatrix + PCMAdjoint);
	// Diagonalize and print out some useful info.
	// This should be done only when debugging...
	/* 
    	std::ofstream matrixOut("PCM_matrix");
    	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(PCMMatrix);
    	if (solver.info() != Eigen::Success)
	       throw std::runtime_error("PCMMatrix diagonalization not successful");	
    	matrixOut << "PCM matrix printout" << std::endl;
 	matrixOut << "Number of Tesserae: " << cavitySize << std::endl;
    	matrixOut << "Largest Eigenvalue: " << solver.eigenvalues()[cavitySize-1] << std::endl;
    	matrixOut << "Lowest Eigenvalue: " << solver.eigenvalues()[0] << std::endl;
    	matrixOut << "Average of Eigenvalues: " << (solver.eigenvalues().sum() / cavitySize)<< std::endl;
    	matrixOut << "List of Eigenvalues:\n" << solver.eigenvalues() << std::endl;
    	matrixOut.close();
	*/
	// Pack into a BlockDiagonalMatrix

	builtIsotropicMatrix = true;
	builtAnisotropicMatrix = false;
}

void IEFSolver::compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge, int irrep) 
{
	if (builtIsotropicMatrix or builtAnisotropicMatrix) 
	{
		charge = - PCMMatrix * potential; // This should look like charge = -PCMMatrix[irrep] * potential; i.e. just using the desired block
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

