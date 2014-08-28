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

#include "IEFSolver.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Cavity.hpp"
#include "IGreensFunction.hpp"
#include "MathUtils.hpp"

void IEFSolver::buildSystemMatrix(const Cavity & cavity)
{
    if (greenInside_->isUniform() && greenOutside_->isUniform()) {
        buildIsotropicMatrix(cavity);
    } else {
        buildAnisotropicMatrix(cavity);
    }
}

void IEFSolver::buildAnisotropicMatrix(const Cavity & cav)
{
    // The total size of the cavity
    int cavitySize = cav.size();
    // The number of irreps in the group
    int nrBlocks = cav.pointGroup().nrIrrep();
    // The size of the irreducible portion of the cavity
    int dimBlock = cav.irreducible_size();

    Eigen::MatrixXd SI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    Eigen::MatrixXd SE = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    Eigen::MatrixXd DI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    Eigen::MatrixXd DE = Eigen::MatrixXd::Zero(cavitySize, cavitySize);

    // This is the very core of PCMSolver
    throw std::runtime_error("Calculation of anisotropic matrix elements not yet implemented.");
/*  greenInside_->compOffDiagonal(cav.elementCenter(), cav.elementNormal(), SI, DI);
    greenInside_->compDiagonal(cav.elementArea(), cav.elementRadius(), SI, DI);
    greenOutside_->compOffDiagonal(cav.elementCenter(), cav.elementNormal(), SE, DE);
    greenOutside_->compDiagonal(cav.elementArea(), cav.elementRadius(), SE, DE);*/
    // Perform symmetry blocking
    // If the group is C1 avoid symmetry blocking, we will just pack the fullPCMMatrix
    // into "block diagonal" when all other manipulations are done.
    if (cav.pointGroup().nrGenerators() != 0) {
        symmetryBlocking(DI, cavitySize, dimBlock, nrBlocks);
        symmetryBlocking(SI, cavitySize, dimBlock, nrBlocks);
        symmetryBlocking(DE, cavitySize, dimBlock, nrBlocks);
        symmetryBlocking(SE, cavitySize, dimBlock, nrBlocks);
    }

    Eigen::MatrixXd a = cav.elementArea().asDiagonal();
    Eigen::MatrixXd aInv = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    aInv = 2 * M_PI * a.inverse();

    fullPCMMatrix = ((aInv - DE) * a * SI + SE * a * (aInv + DI.transpose().eval()));
    fullPCMMatrix = fullPCMMatrix.inverse();
    // 2. Invert T using LU decomposition with full pivoting
    //    This is a rank-revealing LU decomposition, this allows us
    //    to test if T is invertible before attempting to invert it.
    Eigen::FullPivLU<Eigen::MatrixXd> T_LU(fullPCMMatrix);
    if (!(T_LU.isInvertible()))
        throw std::runtime_error("T matrix is not invertible!");
    fullPCMMatrix = T_LU.inverse();
    fullPCMMatrix *= ((aInv - DE) - SE * SI.inverse() * (aInv - DI));
    fullPCMMatrix = fullPCMMatrix * a;
    // 5. Symmetrize K := (K + K+)/2
    if (hermitivitize_) {
        hermitivitize(fullPCMMatrix);
    }
    // Pack into a BlockDiagonalMatrix
    // For the moment just packs into a std::vector<Eigen::MatrixXd> without the syntactic
    // sugar of the BlockDiagonalMatrix class...
    symmetryPacking(blockPCMMatrix, fullPCMMatrix, dimBlock, nrBlocks);
    std::ofstream matrixOut("PCM_matrix");
    matrixOut << "fullPCMMatrix" << std::endl;
    matrixOut << fullPCMMatrix << std::endl;
    for (int i = 0; i < nrBlocks; ++i) {
        matrixOut << "Block number " << i << std::endl;
        matrixOut << blockPCMMatrix[i] << std::endl;
    }

    builtAnisotropicMatrix = true;
    builtIsotropicMatrix = false;
}


void IEFSolver::buildIsotropicMatrix(const Cavity & cav)
{
    // The total size of the cavity
    int cavitySize = cav.size();
    // The number of irreps in the group
    int nrBlocks = cav.pointGroup().nrIrrep();
    // The size of the irreducible portion of the cavity
    int dimBlock = cav.irreducible_size();

    Eigen::MatrixXd SI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    Eigen::MatrixXd DI = Eigen::MatrixXd::Zero(cavitySize, cavitySize);

    // Compute SI and DI on the whole cavity, regardless of symmetry
    for (int i = 0; i < cavitySize; ++i) {
        SI(i, i) = greenInside_->diagonalS(cav.elementArea(i));
        DI(i, i) = greenInside_->diagonalD(cav.elementArea(i), cav.elementRadius(i));
        Eigen::Vector3d source = cav.elementCenter().col(i);
        for (int j = 0; j < cavitySize; ++j) {
            Eigen::Vector3d probe = cav.elementCenter().col(j);
            Eigen::Vector3d probeNormal = cav.elementNormal().col(j);
            probeNormal.normalize();
            if (i != j) {
                SI(i, j) = greenInside_->function(source, probe);
                DI(i, j) = greenInside_->derivative(probeNormal, source, probe);
            }
        }
    }
    // Perform symmetry blocking
    // If the group is C1 avoid symmetry blocking, we will just pack the fullPCMMatrix
    // into "block diagonal" when all other manipulations are done.
    if (cav.pointGroup().nrGenerators() != 0) {
        symmetryBlocking(DI, cavitySize, dimBlock, nrBlocks);
        symmetryBlocking(SI, cavitySize, dimBlock, nrBlocks);
    }

    Eigen::MatrixXd a = cav.elementArea().asDiagonal();
    Eigen::MatrixXd aInv = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    aInv = a.inverse();

    // Tq = -Rv -> q = -(T^-1 * R)v = -Kv
    // T = (2 * M_PI * fact * aInv - DI) * a * SI; R = (2 * M_PI * aInv - DI)
    // fullPCMMatrix = K = T^-1 * R * a
    // 1. Form T
    double epsilon = greenOutside_->epsilon();
    double fact = (epsilon + 1.0)/(epsilon - 1.0);
    fullPCMMatrix = (2 * M_PI * fact * aInv - DI) * a * SI;
    // 2. Invert T using LU decomposition with full pivoting
    //    This is a rank-revealing LU decomposition, this allows us
    //    to test if T is invertible before attempting to invert it.
    Eigen::FullPivLU<Eigen::MatrixXd> T_LU(fullPCMMatrix);
    if (!(T_LU.isInvertible()))
        throw std::runtime_error("T matrix is not invertible!");
    fullPCMMatrix = T_LU.inverse();
    // 3. Multiply T^-1 and R
    fullPCMMatrix *= (2 * M_PI * aInv - DI);
    // 4. Multiply by a
    fullPCMMatrix *= a;
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
    std::ofstream matrixOut("PCM_matrix");
    matrixOut << "fullPCMMatrix" << std::endl;
    matrixOut << fullPCMMatrix << std::endl;
    for (int i = 0; i < nrBlocks; ++i) {
        matrixOut << "Block number " << i << std::endl;
        matrixOut << blockPCMMatrix[i] << std::endl;
    }

    builtIsotropicMatrix = true;
    builtAnisotropicMatrix = false;
}

void IEFSolver::compCharge(const Eigen::VectorXd & potential,
                           Eigen::VectorXd & charge, int irrep)
{
    // The potential and charge vector are of dimension equal to the
    // full dimension of the cavity. We have to select just the part
    // relative to the irrep needed.
    int fullDim = fullPCMMatrix.rows();
    int nrBlocks = blockPCMMatrix.size();
    int irrDim = int(fullDim/nrBlocks);
    if (builtIsotropicMatrix or builtAnisotropicMatrix) {
        charge.segment(irrep*irrDim,
                       irrDim)= - blockPCMMatrix[irrep] * potential.segment(irrep*irrDim, irrDim);
//		charge = - fullPCMMatrix * potential;
    } else {
        throw std::runtime_error("PCM matrix not initialized!");
    }
}

std::ostream & IEFSolver::printSolver(std::ostream & os)
{
    std::string type;
    if (builtAnisotropicMatrix) {
        type = "IEFPCM, anisotropic";
    } else {
        type = "IEFPCM, isotropic";
    }
    os << "Solver Type: " << type << std::endl;
    if (hermitivitize_) {
        os << "PCM matrix hermitivitized";
    } else {
        os << "PCM matrix NOT hermitivitized (matches old DALTON)";
    }
    return os;
}

