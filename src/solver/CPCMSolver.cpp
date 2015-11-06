/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#include "CPCMSolver.hpp"

#include <fstream>
#include <iostream>
#include <string>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include "Cavity.hpp"
#include "Element.hpp"
#include "IGreensFunction.hpp"
#include "MathUtils.hpp"

void CPCMSolver::buildSystemMatrix_impl(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o)
{
    if (!isotropic_) PCMSOLVER_ERROR("C-PCM is defined only for isotropic environments!");

    // The total size of the cavity
    size_t cavitySize = cavity.size();
    // The number of irreps in the group
    int nrBlocks = cavity.pointGroup().nrIrrep();
    // The size of the irreducible portion of the cavity
    int dimBlock = cavity.irreducible_size();

    // Compute SI on the whole cavity, regardless of symmetry
    Eigen::MatrixXd SI = gf_i.singleLayer(cavity.elements());

    // Perform symmetry blocking only for the SI matrix as the DI matrix is not used.
    // If the group is C1 avoid symmetry blocking, we will just pack the fullPCMMatrix
    // into "block diagonal" when all other manipulations are done.
    if (cavity.pointGroup().nrGenerators() != 0) {
        symmetryBlocking(SI, cavitySize, dimBlock, nrBlocks);
    }

    double epsilon = profiles::epsilon(gf_o.permittivity());
    double fact = (epsilon - 1.0)/(epsilon + correction_);
    // Invert SI  using LU decomposition with full pivoting
    // This is a rank-revealing LU decomposition, this allows us
    // to test if SI is invertible before attempting to invert it.
    Eigen::FullPivLU<Eigen::MatrixXd> SI_LU(SI);
    if (!(SI_LU.isInvertible()))
        PCMSOLVER_ERROR("SI matrix is not invertible!");
    fullPCMMatrix_ = fact * SI_LU.inverse();
    // 5. Symmetrize K := (K + K+)/2
    if (hermitivitize_) {
        hermitivitize(fullPCMMatrix_);
    }
    symmetryPacking(blockPCMMatrix_, fullPCMMatrix_, dimBlock, nrBlocks);

    built_ = true;
}

Eigen::VectorXd CPCMSolver::computeCharge_impl(const Eigen::VectorXd & potential, int irrep) const
{
    // The potential and charge vector are of dimension equal to the
    // full dimension of the cavity. We have to select just the part
    // relative to the irrep needed.
    int fullDim = fullPCMMatrix_.rows();
    Eigen::VectorXd charge = Eigen::VectorXd::Zero(fullDim);
    int nrBlocks = blockPCMMatrix_.size();
    int irrDim = fullDim/nrBlocks;
    charge.segment(irrep*irrDim, irrDim) =
        - blockPCMMatrix_[irrep] * potential.segment(irrep*irrDim, irrDim);
    return charge;
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

