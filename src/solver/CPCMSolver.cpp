#include "CPCMSolver.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "Cavity.hpp"
#include "GreensFunction.hpp"
#include "MathUtils.hpp"

void CPCMSolver::buildSystemMatrix(const Cavity & cavity)
{
    if (greenInside_->isUniform() && greenOutside_->isUniform()) {
        buildIsotropicMatrix(cavity);
    } else {
        throw std::runtime_error("C-PCM is defined only for isotropic environments!");
    }
}


void CPCMSolver::buildIsotropicMatrix(const Cavity & cav)
{
    // The total size of the cavity
    int cavitySize = cav.size();
    // The number of irreps in the group
    int nrBlocks = cav.pointGroup().nrIrrep();
    // The size of the irreducible portion of the cavity
    int dimBlock = cav.irreducible_size();

    Eigen::MatrixXd SI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    Eigen::MatrixXd DI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);

    // This is the very core of PCMSolver
    greenInside_->compOffDiagonal(cav.elementCenter(), cav.elementNormal(), SI, DI);
    greenInside_->compDiagonal(cav.elementArea(), cav.elementRadius(), SI, DI);
    // Perform symmetry blocking only for the SI matrix as the DI matrix is not used.
    // If the group is C1 avoid symmetry blocking, we will just pack the fullPCMMatrix
    // into "block diagonal" when all other manipulations are done.
    if (cav.pointGroup().nrGenerators() != 0) {
        symmetryBlocking(SI, cavitySize, dimBlock, nrBlocks);
    }

    double epsilon = greenOutside_->dielectricConstant();
    double fact = (epsilon - 1.0)/(epsilon + correction_);
    // Invert SI  using LU decomposition with full pivoting
    // This is a rank-revealing LU decomposition, this allows us
    // to test if SI is invertible before attempting to invert it.
    Eigen::FullPivLU<Eigen::MatrixXd> SI_LU(SI);
    if (!(SI_LU.isInvertible()))
        throw std::runtime_error("SI matrix is not invertible!");
    fullPCMMatrix = fact * SI_LU.inverse();
    // 5. Symmetrize K := (K + K+)/2
    if (hermitivitize_) {
        hermitivitize(fullPCMMatrix);
    }
    // Diagonalize and print out some useful info.
    // This should be done only when debugging...
    /*
    	std::ofstream matrixOut("PCM_matrix");
    	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(fullPCMMatrix);
    	if (solver.info() != Eigen::Success)
           throw std::runtime_error("fullPCMMatrix diagonalization not successful");
    	matrixOut << "PCM matrix printout" << std::endl;
    matrixOut << "Number of Tesserae: " << cavitySize << std::endl;
    	matrixOut << "Largest Eigenvalue: " << solver.eigenvalues()[cavitySize-1] << std::endl;
    	matrixOut << "Lowest Eigenvalue: " << solver.eigenvalues()[0] << std::endl;
    	matrixOut << "Average of Eigenvalues: " << (solver.eigenvalues().sum() / cavitySize)<< std::endl;
    	matrixOut << "List of Eigenvalues:\n" << solver.eigenvalues() << std::endl;
    	matrixOut.close();
    */
    // Pack into a BlockDiagonalMatrix
    // For the moment just packs into a std::vector<Eigen::MatrixXd> without the syntactic
    // sugar of the BlockDiagonalMatrix class...
    symmetryPacking(blockPCMMatrix, fullPCMMatrix, dimBlock, nrBlocks);

    builtIsotropicMatrix = true;
    builtAnisotropicMatrix = false;
}

void CPCMSolver::compCharge(const Eigen::VectorXd & potential,
                            Eigen::VectorXd & charge, int irrep)
{
    // The potential and charge vector are of dimension equal to the
    // full dimension of the cavity. We have to select just the part
    // relative to the irrep needed.
    int fullDim = fullPCMMatrix.rows();
    int nrBlocks = blockPCMMatrix.size();
    int irrDim = int(fullDim/nrBlocks);
    if (builtIsotropicMatrix) {
        charge.segment(irrep*irrDim,
                       irrDim)= - blockPCMMatrix[irrep] * potential.segment(irrep*irrDim, irrDim);
//		charge = - fullPCMMatrix * potential;
    } else {
        throw std::runtime_error("PCM matrix not initialized!");
    }
}

std::ostream & CPCMSolver::printSolver(std::ostream & os)
{
    os << "Solver Type: C-PCM" << std::endl;
    if (hermitivitize_) {
        os << "PCM matrix hermitivitized";
    } else {
        os << "PCM matrix NOT hermitivitized (matches old DALTON)";
    }
    return os;
}

