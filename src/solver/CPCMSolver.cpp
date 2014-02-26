#include "CPCMSolver.hpp"

#include <fstream>
#include <ostream>
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

void CPCMSolver::buildSystemMatrix(const Cavity & cavity) 
{
	if (greenInside_->isUniform() && greenOutside_->isUniform()) 
	{
		buildIsotropicMatrix(cavity);
	} 
	else 
	{
	        throw std::runtime_error("C-PCM is defined only for isotropic environments!");
	}
}


void CPCMSolver::buildIsotropicMatrix(const Cavity & cav)
{
	bool symmetry = cav.pointGroup().groupInteger();
	double epsilon = greenOutside_->getDielectricConstant();
    	int cavitySize = cav.size();
    	Eigen::MatrixXd SI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
	Eigen::MatrixXd DI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    	
	// This is the very core of PCMSolver
    	greenInside_->compOffDiagonal(cav.getElementCenter(), cav.getElementNormal(), SI, DI);
    	greenInside_->compDiagonal(cav.getElementArea(), cav.getElementRadius(), SI, DI);
	// Perform symmetry blocking only for the SI matrix as the DI matrix is not used.
	// If the group is C1 avoid symmetry blocking, we will just pack the PCMMatrix
	// into "block diagonal" when all other manipulations are done.
	if (symmetry)
	{
		symmetryBlocking(SI, cav);
	}

	Eigen::MatrixXd a = cav.getElementArea().asDiagonal();    
	
	double fact = (epsilon - 1.0)/(epsilon + correction_);
	// Invert SI  using LU decomposition with full pivoting
	// This is a rank-revealing LU decomposition, this allows us
	// to test if SI is invertible before attempting to invert it.
	Eigen::FullPivLU<Eigen::MatrixXd> SI_LU(SI);
	if (!(SI_LU.isInvertible()))
		throw std::runtime_error("SI matrix is not invertible!");
    	PCMMatrix = fact * SI_LU.inverse();
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

void CPCMSolver::compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge, int irrep) 
{
	if (builtIsotropicMatrix) 
	{
		charge = - PCMMatrix * potential;
	} 
	else 
	{
		throw std::runtime_error("PCM matrix not initialized!");
	}
}
    
std::ostream & CPCMSolver::printSolver(std::ostream & os) 
{
	os << "Solver Type: C-PCM";
	return os;
}

