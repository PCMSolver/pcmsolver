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

#include "IEFSolver.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include "cavity/Cavity.hpp"
#include "cavity/Element.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"
#include "SolverImpl.hpp"

void IEFSolver::buildSystemMatrix_impl(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o)
{
    isotropic_ = (gf_i.uniform() && gf_o.uniform());
    isotropic_ ? buildIsotropicMatrix(cavity, gf_i, gf_o) : buildAnisotropicMatrix(cavity, gf_i, gf_o);
}

void IEFSolver::buildAnisotropicMatrix(const Cavity & cav, const IGreensFunction & gf_i, const IGreensFunction & gf_o)
{
    fullPCMMatrix_ = anisotropicIEFMatrix(cav, gf_i, gf_o);
    // Symmetrize K := (K + K+)/2
    if (hermitivitize_) {
        hermitivitize(fullPCMMatrix_);
    }
    // Pack into a block diagonal matrix
    // The number of irreps in the group
    int nrBlocks = cav.pointGroup().nrIrrep();
    // The size of the irreducible portion of the cavity
    int dimBlock = cav.irreducible_size();
    // For the moment just packs into a std::vector<Eigen::MatrixXd>
    symmetryPacking(blockPCMMatrix_, fullPCMMatrix_, dimBlock, nrBlocks);

    built_ = true;
}

void IEFSolver::buildIsotropicMatrix(const Cavity & cav, const IGreensFunction & gf_i, const IGreensFunction & gf_o)
{
    fullPCMMatrix_ = isotropicIEFMatrix(cav, gf_i, profiles::epsilon(gf_o.permittivity()));
    // Symmetrize K := (K + K+)/2
    if (hermitivitize_) {
        hermitivitize(fullPCMMatrix_);
    }
    // Pack into a block diagonal matrix
    // The number of irreps in the group
    int nrBlocks = cav.pointGroup().nrIrrep();
    // The size of the irreducible portion of the cavity
    int dimBlock = cav.irreducible_size();
    // For the moment just packs into a std::vector<Eigen::MatrixXd>
    symmetryPacking(blockPCMMatrix_, fullPCMMatrix_, dimBlock, nrBlocks);

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
    int irrDim = fullDim/nrBlocks;
    charge.segment(irrep*irrDim, irrDim) =
        - blockPCMMatrix_[irrep] * potential.segment(irrep*irrDim, irrDim);
    return charge;
}

std::ostream & IEFSolver::printSolver(std::ostream & os)
{
    std::string type;
    if (isotropic_) {
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

