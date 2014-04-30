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

#include "WEMSolver.hpp"

#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include "Config.hpp"

#include "EigenPimpl.hpp"

extern "C"
{
#include "vector3.h"
#include "sparse2.h"
#include "intvector.h"
#include "basis.h"
#include "WEM.h"
#include "read_points.h"
#include "vector2.h"
#include "interpolate.h"
#include "topology.h"
#include "kern.h"
#include "compression.h"
#include "postproc.h"
#include "WEMRHS.h"
#include "WEMPCG.h"
#include "WEMPGMRES.h"
#include "dwt.h"
#include "cubature.h"
#include "gauss_square.h"
#include "constants.h"
}

#include "Cavity.hpp"
#include "IGreensFunction.hpp"
#include "WaveletCavity.hpp"

void WEMSolver::initWEMMembers()
{
    pointList = NULL;
    nodeList = NULL;
    elementList = NULL;
    T_ = NULL;
    systemMatricesInitialized_ = false;
    threshold = 1e-10;
    quadratureLevel_ = 1;
    nQuadPoints = 0;
}

WEMSolver::~WEMSolver()
{
    if (nodeList != NULL)    free(nodeList);
    if (elementList != NULL) free_patchlist(&elementList,nFunctions);
    if (pointList != NULL)   free_points(&pointList, nPatches, nLevels);
    if (systemMatricesInitialized_) {
        free_sparse2(&S_i_);
        if(integralEquation == Full) free_sparse2(&S_e_);
    }
}

void WEMSolver::uploadCavity(const WaveletCavity & cavity)
{
    nPatches = cavity.getNPatches();
    nLevels = cavity.getNLevels();
    int n = (1<<nLevels);
    nFunctions = nPatches * n * n;
    alloc_points(&pointList, nPatches, nLevels);
    int kk = 0;
    // Ask Helmut about index switch
    for (size_t i = 0; i < nPatches; ++i) {
        for (int j = 0; j <= n; ++j) {
            for (int k = 0; k <= n; ++k) {
                Eigen::Vector3d p = cavity.getNodePoint(kk);
                pointList[i][k][j] = vector3_make(p(0), p(1), p(2));
                kk++;
            }
        }
    }
}

void WEMSolver::buildSystemMatrix(const Cavity & cavity)
{
    // Down-cast const Cavity & to const WaveletCavity &
    // This is messy. The wavelet classes are not really Object-Oriented...
    // In principle, I should be able to use Cavity directly, without
    // traversing the inheritance hierarchy!
    try {
        const WaveletCavity & waveletCavity = dynamic_cast<const WaveletCavity&>(cavity);
        uploadCavity(waveletCavity);
        initInterpolation();
        constructWavelets();
        constructSystemMatrix();
    } catch (const std::bad_cast & e) {
        throw std::runtime_error(e.what() +
                                 std::string(" Wavelet type cavity needed for wavelet solver."));
    }
}

void WEMSolver::constructSystemMatrix()
{
    constructSi();
    if(integralEquation == Full) {
        constructSe();
    }
}

void WEMSolver::compCharge(const Eigen::VectorXd & potential,
                           Eigen::VectorXd & charge, int irrep)
{
    switch (integralEquation) {
    case FirstKind:
        solveFirstKind(potential, charge);
        break;
    case SecondKind:
        solveSecondKind(potential, charge);
        break;
    case Full:
        solveFull(potential, charge);
        break;
    default:
        throw std::runtime_error("Invalid case");
    }
    charge *= -1.0;
    //	charge /= -ToAngstrom; //WARNING  WARNING  WARNING
}
