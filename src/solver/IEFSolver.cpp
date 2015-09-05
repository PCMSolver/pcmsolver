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
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#include "IEFSolver.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include "Cavity.hpp"
#include "Element.hpp"
#include "IGreensFunction.hpp"
#include "MathUtils.hpp"

void IEFSolver::buildSystemMatrix_impl(const Cavity & cavity)
{
    if (greenInside_->uniform() && greenOutside_->uniform()) {
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

    // Compute SI, DI and SE, DE on the whole cavity, regardless of symmetry
    Eigen::MatrixXd SI = greenInside_->singleLayer(cav.elements());
    Eigen::MatrixXd DI = greenInside_->doubleLayer(cav.elements());
    Eigen::MatrixXd SE = greenOutside_->singleLayer(cav.elements());
    Eigen::MatrixXd DE = greenOutside_->doubleLayer(cav.elements());

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
    Eigen::MatrixXd aInv = a.inverse();

    // 1. Form T
    fullPCMMatrix_ = ((2 * M_PI * aInv - DE) * a * SI + SE * a * (2 * M_PI * aInv + DI.transpose().eval()));
    // 2. Invert T using LU decomposition with full pivoting
    //    This is a rank-revealing LU decomposition, this allows us
    //    to test if T is invertible before attempting to invert it.
    Eigen::FullPivLU<Eigen::MatrixXd> T_LU(fullPCMMatrix_);
    if (!(T_LU.isInvertible())) PCMSOLVER_ERROR("T matrix is not invertible!");
    fullPCMMatrix_ = T_LU.inverse();
    Eigen::FullPivLU<Eigen::MatrixXd> SI_LU(SI);
    if (!(SI_LU.isInvertible())) PCMSOLVER_ERROR("SI matrix is not invertible!");
    fullPCMMatrix_ *= ((2 * M_PI * aInv - DE) - SE * SI_LU.inverse() * (2 * M_PI * aInv - DI));
    fullPCMMatrix_ *= a;
    // 5. Symmetrize K := (K + K+)/2
    if (hermitivitize_) {
        hermitivitize(fullPCMMatrix_);
    }
    // Pack into a block diagonal matrix
    // For the moment just packs into a std::vector<Eigen::MatrixXd>
    symmetryPacking(blockPCMMatrix_, fullPCMMatrix_, dimBlock, nrBlocks);
    std::ofstream matrixOut("PCM_matrix");
    matrixOut << "fullPCMMatrix" << std::endl;
    matrixOut << fullPCMMatrix_ << std::endl;
    for (int i = 0; i < nrBlocks; ++i) {
        matrixOut << "Block number " << i << std::endl;
        matrixOut << blockPCMMatrix_[i] << std::endl;
    }

    built_ = true;
}

void IEFSolver::buildIsotropicMatrix(const Cavity & cav)
{
    // The total size of the cavity
    int cavitySize = cav.size();
    // The number of irreps in the group
    int nrBlocks = cav.pointGroup().nrIrrep();
    // The size of the irreducible portion of the cavity
    int dimBlock = cav.irreducible_size();

    // Compute SI and DI on the whole cavity, regardless of symmetry
    Eigen::MatrixXd SI = greenInside_->singleLayer(cav.elements());
    Eigen::MatrixXd DI = greenInside_->doubleLayer(cav.elements());

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
    // fullPCMMatrix_ = K = T^-1 * R * a
    // 1. Form T
    double epsilon = profiles::epsilon(greenOutside_->permittivity());
    double fact = (epsilon + 1.0)/(epsilon - 1.0);
    fullPCMMatrix_ = (2 * M_PI * fact * aInv - DI) * a * SI;
    // 2. Invert T using LU decomposition with full pivoting
    //    This is a rank-revealing LU decomposition, this allows us
    //    to test if T is invertible before attempting to invert it.
    Eigen::FullPivLU<Eigen::MatrixXd> T_LU(fullPCMMatrix_);
    if (!(T_LU.isInvertible()))
        PCMSOLVER_ERROR("T matrix is not invertible!");
    fullPCMMatrix_ = T_LU.inverse();
    // 3. Multiply T^-1 and R
    fullPCMMatrix_ *= (2 * M_PI * aInv - DI);
    // 4. Multiply by a
    fullPCMMatrix_ *= a;
    // 5. Symmetrize K := (K + K+)/2
    if (hermitivitize_) {
        hermitivitize(fullPCMMatrix_);
    }
    // Pack into a block diagonal matrix
    // For the moment just packs into a std::vector<Eigen::MatrixXd>
    symmetryPacking(blockPCMMatrix_, fullPCMMatrix_, dimBlock, nrBlocks);
    std::ofstream matrixOut("PCM_matrix");
    matrixOut << "fullPCMMatrix" << std::endl;
    matrixOut << fullPCMMatrix_ << std::endl;
    for (int i = 0; i < nrBlocks; ++i) {
        matrixOut << "Block number " << i << std::endl;
        matrixOut << blockPCMMatrix_[i] << std::endl;
    }

    built_ = true;
}

Eigen::VectorXd IEFSolver::computeCharge_impl(const Eigen::VectorXd & potential, int irrep) const
{
    // The potential and charge vector are of dimension equal to the
    // full dimension of the cavity. We have to select just the part
    // relative to the irrep needed.
    int fullDim = fullPCMMatrix_.rows();
    Eigen::VectorXd charge = Eigen::VectorXd::Zero(fullDim);
    int nrBlocks = blockPCMMatrix_.size();
    int irrDim = int(fullDim/nrBlocks);
    charge.segment(irrep*irrDim, irrDim) =
        - blockPCMMatrix_[irrep] * potential.segment(irrep*irrDim, irrDim);
    return charge;
}

std::ostream & IEFSolver::printSolver(std::ostream & os)
{
    std::string type;
    if (greenInside_->uniform() && greenOutside_->uniform()) {
        type = "IEFPCM, isotropic";
    } else {
        type = "IEFPCM, anisotropic";
    }
    os << "Solver Type: " << type << std::endl;
    if (hermitivitize_) {
        os << "PCM matrix hermitivitized";
    } else {
        os << "PCM matrix NOT hermitivitized (matches old DALTON)";
    }
    return os;
}

