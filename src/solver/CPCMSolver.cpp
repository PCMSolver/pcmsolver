/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *     
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify       
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *                                                                          
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *                                                                          
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#include "CPCMSolver.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "Cavity.hpp"
#include "Element.hpp"
#include "IGreensFunction.hpp"
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

    // Compute SI on the whole cavity, regardless of symmetry
    for (int i = 0; i < cavitySize; ++i) {
        SI(i, i) = greenInside_->diagonalS(cav.elements(i));
        Eigen::Vector3d source = cav.elementCenter().col(i);
        for (int j = 0; j < cavitySize; ++j) {
            Eigen::Vector3d probe = cav.elementCenter().col(j);
            Eigen::Vector3d probeNormal = cav.elementNormal().col(j);
            probeNormal.normalize();
            if (i != j) {
                SI(i, j) = greenInside_->function(source, probe);
            }
        }
    }
    // Perform symmetry blocking only for the SI matrix as the DI matrix is not used.
    // If the group is C1 avoid symmetry blocking, we will just pack the fullPCMMatrix
    // into "block diagonal" when all other manipulations are done.
    if (cav.pointGroup().nrGenerators() != 0) {
        symmetryBlocking(SI, cavitySize, dimBlock, nrBlocks);
    }

    double epsilon = greenOutside_->epsilon();
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

