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

#include "WEMSolver.hpp"

#include <iostream>
#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

#include "WEM.hpp"
#include "WEMPGMRES.hpp"
#include "WEMPCG.hpp"
#include "WEMRHS.hpp"
#include "Energy.hpp"

#include "Cavity.hpp"
#include "IGreensFunction.hpp"
#include "WaveletCavity.hpp"

static IGreensFunction *gf;

static double SingleLayer (Vector3 x, Vector3 y)
{
    Eigen::Vector3d vx(x.x, x.y, x.z);
    Eigen::Vector3d vy(y.x, y.y, y.z);
    Eigen::Vector3d foo = Eigen::Vector3d::Zero();
    double value = gf->function(vx, vy);
    return value;
}

static double DoubleLayer (Vector3 x, Vector3 y, Vector3 n_y)
{
    Eigen::Vector3d vx(x.x, x.y, x.z);
    Eigen::Vector3d vy(y.x, y.y, y.z);
    Eigen::Vector3d vn_y(n_y.x, n_y.y, n_y.z);
    Eigen::Vector3d foo = Eigen::Vector3d::Zero();
    double value = gf->derivative(vn_y, vx, vy);
    return value;
}

static double SLInt(Vector3 x, Vector3 y)
{
    return(1.0l/sqrt((x.x-y.x)*(x.x-y.x)+(x.y-y.y)*(x.y-y.y)+(x.z-y.z)*(x.z-y.z)));
}

static double SLExt(Vector3 x, Vector3 y)
{
    double r = sqrt((x.x-y.x)*(x.x-y.x)+(x.y-y.y)*(x.y-y.y)+(x.z-y.z)*(x.z-y.z));
    return 1.0l/(r*78.39l);
}

static double DLUni(Vector3 x, Vector3 y, Vector3 n_y)
{
    Vector3	c;
    double r;
    Eigen::Vector3d grad, dir;
    c.x = x.x-y.x;
    c.y = x.y-y.y;
    c.z = x.z-y.z;
    r = sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
    grad(0) = c.x/(r*r*r);
    grad(1) = c.y/(r*r*r);
    grad(2) = c.z/(r*r*r);
    dir << n_y.x, n_y.y, n_y.z;
    return (c.x*n_y.x+c.y*n_y.y+c.z*n_y.z)/(r*r*r);
}

void WEMSolver::initWEMMembers()
{
    pointList = NULL;
//    nodeList = NULL;
//    elementList = NULL;
//    T_ = NULL;
    systemMatricesInitialized_ = false;
    threshold = 1e-10;
    //af->quadratureLevel_ = 1; // set in constructor of AnsatzFunction
    //af->nQuadPoints = 0; //??? what is this for?
    af->elementTree.element = NULL;
    af->elementTree.nop = 0;
    af->waveletList.W = NULL;
    af->waveletList.sizeWaveletList = 0;
}

WEMSolver::~WEMSolver()
{
    // TODO insert delete of af - ansatzFunction
    delete(af);
    //if (nodeList != NULL)    free(nodeList);
    //if (elementList != NULL) free_patchlist(&elementList,nFunctions);
    if (pointList != NULL)   freePoints(&pointList);
    if (systemMatricesInitialized_) {
        freeSparse2(&S_i_);
        if(integralEquation == Full) freeSparse2(&S_e_);
    }
    //if(elementTree != NULL) free_elementlist(&elementTree, nPatches,nLevels);
    //if(waveletList != NULL) free_waveletlist(&waveletList, nPatches,nLevels);
    //if(T_ != NULL)          free_interpolate(&T_,nPatches,nLevels);
}

void WEMSolver::uploadCavity(const WaveletCavity & cavity)
{
    af->nPatches = cavity.getNPatches();
    af->nLevels = cavity.getNLevels();
    int n = (1<<af->nLevels);
    af->nFunctions = af->nPatches * n * n;
    allocPoints(&pointList, af->nPatches, af->nLevels);
    int kk = 0;
    // index switch for "faster access"
    for (size_t i = 0; i < af->nPatches; ++i) {
        for (int j = 0; j <= n; ++j) {
            for (int k = 0; k <= n; ++k) {
                Eigen::Vector3d p = cavity.getNodePoint(kk);
                pointList[i][k][j] = Vector3(p(0), p(1), p(2));
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

void WEMSolver::initInterpolation(){
    af->interCoeff = new Interpolation(pointList, interpolationGrade, interpolationType, af->nLevels, af->nPatches);
    af->nNodes = af->genNet(pointList);
}

void WEMSolver::constructWavelets(){
    //af->generateElementList(); // already done in genNet
    af->generateWaveletList();
    af->setQuadratureLevel();
    af->simplifyWaveletList();
    //af->completeElementList();
}

void WEMSolver::constructSystemMatrix()
{
    constructSi();
    if(integralEquation == Full) {
        constructSe();
    }
}

void WEMSolver::constructSi(){
    double factor = 0.0;
    double epsilon = 0;
    switch (integralEquation) {
        case FirstKind:
            epsilon = greenOutside_->epsilon();
            factor = - 2 * M_PI * (epsilon + 1) / (epsilon - 1);
            break;
        case SecondKind:
            throw std::runtime_error("Second Kind not yet implemented"); //careful to the double layer sign when implementing it....
            break;
        case Full:
            factor = 2 * M_PI;
            break;
        default:
            throw std::runtime_error("Unknown integral equation type.");
    }
    gf = greenInside_;
    apriori1_ = af->compression(&S_i_);
    WEM(af, &S_i_, SingleLayer, DoubleLayer, factor);
    aposteriori1_ = af->postProc(&S_i_);
}

void WEMSolver::constructSe(){
    gf = greenOutside_; // sets the global pointer to pass GF to C code
    apriori2_ = af->compression(&S_e_);
    WEM(af, &S_e_, SingleLayer, DoubleLayer, -2*M_PI);
    aposteriori2_ = af->postProc(&S_e_);
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

void WEMSolver::solveFirstKind(const Eigen::VectorXd & potential,
                               Eigen::VectorXd & charge)
{
    double *rhs;
    double *u = (double*) calloc(af->nFunctions, sizeof(double));
    double * pot = const_cast<double *>(potential.data());
    double * chg = charge.data();
    double epsilon = greenOutside_->epsilon();
    WEMRHS2M(&rhs, pot, af);
    int iter = WEMPGMRES2(&S_i_, rhs, u, threshold, af);
    af->tdwt(u);
    af->dwt(u);
    for (size_t i = 0; i < af->nFunctions; ++i) {
        rhs[i] += 4 * M_PI * u[i] / (epsilon - 1);
    }
    memset(u, 0, af->nFunctions * sizeof(double));
    iter = WEMPCG(&S_i_, rhs, u, threshold, af);
    af->tdwt(u);
    energy_ext(u, pot, af);
    charge_ext(u, chg, af);
    free(rhs);
    free(u);
}

void WEMSolver::solveSecondKind(const Eigen::VectorXd & potential,
                                Eigen::VectorXd & charge)
{
    throw std::runtime_error("Second Kind not yet implemented"); //careful to the double layer sign when implementing it....
}

void WEMSolver::solveFull(const Eigen::VectorXd & potential,
                          Eigen::VectorXd & charge)
{
    double *rhs;
    double *u = (double*) calloc(af->nFunctions, sizeof(double));
    double *v = (double*) calloc(af->nFunctions, sizeof(double));
    //next line is just a quick fix to avoid problems with const but i do not like it...
    double * pot = const_cast<double *>(potential.data());
    WEMRHS2M(&rhs,pot, af);
    int iters = WEMPCG(&S_i_, rhs, u, threshold, af);
    memset(rhs, 0, af->nFunctions*sizeof(double));
    for(unsigned int i = 0; i < af->nFunctions; i++) {
        for(unsigned int j = 0; j < S_e_.row_number[i]; j++) {
            rhs[i] += S_e_.value1[i][j] * u[S_e_.index[i][j]];
        }
    }
    iters = WEMPGMRES3(&S_i_, &S_e_, rhs, v, threshold, af);
    for(unsigned int i = 0; i < af->nFunctions; i++) {
        u[i] -= 4*M_PI*v[i];
    }
    af->tdwt(u);
    energy_ext(u, pot, af);
    charge_ext(u, charge.data(), af);
    free(rhs);
    free(u);
    free(v);
}

std::ostream & WEMSolver::printSolver(std::ostream & os)
{
    os << "Solver Type: Wavelet, piecewise constant or linear functions";
    return os;
}
