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

#include "PWLSolver.hpp"

#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

extern "C"
{
#include "vector3.h"
#include "sparse.h"
#include "sparse2.h"
#include "intvector_pwl.h"
#include "basis_pwl.h"
#include "WEM_pwl.h"
#include "read_points.h"
#include "vector2.h"
#include "interpolate_pwl.h"
#include "topology_pwl.h"
#include "kern.h"
#include "compression_pwl.h"
#include "postproc_pwl.h"
#include "WEMRHS_pwl.h"
#include "WEMPCG_pwl.h"
#include "WEMPGMRES_pwl.h"
#include "dwt_pwl.h"
#include "cubature.h"
#include "gauss_square.h"
#include "constants.h"
#include "precond_pwl.h"
#include "energy_pwl.h"
}

#include "Cavity.hpp"
#include "IGreensFunction.hpp"
#include "WaveletCavity.hpp"

static IGreensFunction * gf;

static double SLInt(vector3 x, vector3 y)
{
    return(1.0l/sqrt((x.x-y.x)*(x.x-y.x)+(x.y-y.y)*(x.y-y.y)+(x.z-y.z)*(x.z-y.z)));
}


static double SLExt(vector3 x, vector3 y)
{
    double r = sqrt((x.x-y.x)*(x.x-y.x)+(x.y-y.y)*(x.y-y.y)+(x.z-y.z)*(x.z-y.z));
    return 1.0l/(r*78.39l);
}


static double DLUni(vector3 x, vector3 y, vector3 n_y)
{
    vector3	c;
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

static double SingleLayer (vector3 x, vector3 y)
{
    Eigen::Vector3d vx(x.x, x.y, x.z);
    Eigen::Vector3d vy(y.x, y.y, y.z);
    double value = gf->function(vx, vy);
    return value;
}

static double DoubleLayer (vector3 x, vector3 y, vector3 n_y)
{
    Eigen::Vector3d vx(x.x, x.y, x.z);
    Eigen::Vector3d vy(y.x, y.y, y.z);
    Eigen::Vector3d vn_y(n_y.x, n_y.y, n_y.z);
    double value = gf->derivative(vn_y, vx, vy);
    return value;
}

void PWLSolver::initPointers()
{
    elementTree = NULL;
    waveletList = NULL;
}

/*PWLSolver::PWLSolver(Section solver) : WEMSolver(solver)
  {
	initPointers();
	setSolverType("Linear");
}*/

PWLSolver::~PWLSolver()
{
    if(elementTree != NULL) free_elementlist_pwl(&elementTree, nPatches, nLevels);
    if(waveletList != NULL) free_waveletlist_pwl(&waveletList, nNodes);
    if(T_ != NULL)          free_interpolate_pwl(&T_,nPatches,nLevels);
}

void PWLSolver::initInterpolation()
{
    init_interpolate_pwl(&T_, pointList, nPatches, nLevels);
    nNodes = gennet_pwl(&nodeList, &elementList, pointList, nPatches, nLevels);
FILE* debugFile = fopen("debug.out","w");
fclose(debugFile);
debugFile = fopen("debug.out","a");	
fprintf(debugFile,">>> PPOINTLIST\n");
for(unsigned int m = 0; m<nNodes; ++m){
	fprintf(debugFile,"%lf %lf %lf\n",(nodeList)[m].x, (nodeList)[m].y, (nodeList)[m].z);
}
fprintf(debugFile,"<<< PPOINTLIST\n");
fprintf(debugFile,">>> PELEMENTLIST\n");
unsigned int nf = nPatches*(1<<nLevels)*(1<<nLevels);	        /* Anzahl der Patches */
for(unsigned int m = 0; m<nf; ++m){
	fprintf(debugFile,"%d %d %d %d\n",(elementList)[m][0], (elementList)[m][1], (elementList)[m][2], (elementList)[m][3]);
}
fprintf(debugFile,"<<< PELEMENTLIST\n");
fclose(debugFile);
}

void PWLSolver::constructWavelets()
{
    generate_elementlist_pwl(&elementTree, nodeList, elementList, nPatches, nLevels);
    generate_waveletlist_pwl(&waveletList, elementTree, nPatches, nLevels, nNodes);
    set_quadrature_level_pwl(waveletList, elementTree, nPatches, nLevels, nNodes);
    simplify_waveletlist_pwl(waveletList, elementTree, nPatches, nLevels, nNodes);
    complete_elementlist_pwl(waveletList, elementTree, nPatches, nLevels, nNodes);
unsigned int	N = 1 << nLevels;	/* N*N Elemente pro Patch auf dem Level M     */
unsigned int	ne;		/* Anzahl der Elemente                        */
ne = nPatches*(4*N*N-1)/3;			/* Anzahl der Elemente */
FILE* debugFile = fopen("debug.out","a");
fprintf(debugFile,">>> HIERARCHICAL_ELEMENT_TREE\n");
for(unsigned int m = 0; m<ne; ++m){
    fprintf(debugFile,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf\n", elementTree[m].patch, elementTree[m].level, elementTree[m].index_s, elementTree[m].index_t, elementTree[m].midpoint.x, elementTree[m].midpoint.y, elementTree[m].midpoint.z, elementTree[m].radius, nodeList[elementTree[m].vertex[0]].x, nodeList[elementTree[m].vertex[0]].y, nodeList[elementTree[m].vertex[0]].z);
    for(unsigned int i1 = 0; i1< elementTree[m].wavelet_number;++i1)
        fprintf(debugFile,"%d ", elementTree[m].wavelet[i1]);
    fprintf(debugFile,"\n");
}
fprintf(debugFile,"<<< HIERARCHICAL_ELEMENT_TREE\n");
fclose(debugFile);
}

void PWLSolver::constructSi()
{
    double factor = 0;
    double epsilon = 0;
    switch (integralEquation) {
    case FirstKind:
    case SecondKind:
        epsilon = greenOutside_->epsilon();
        factor = - 2 * M_PI * (epsilon + 1) / (epsilon - 1);
        break;
    case Full:
        factor = 2 * M_PI;
        break;
    default:
        throw std::runtime_error("Unknown integral equation type.");

    }
    gf = greenInside_;
    apriori1_ = compression_pwl(&S_i_, waveletList, elementTree, nPatches, nLevels,
                                nNodes);
FILE* debugFile = fopen("debug.out", "a");
fprintf(debugFile,">>> SYSTEMMATRIX AC\n");
fprintf(debugFile, "%d %d\n", S_i_.m, S_i_.n);
for(unsigned int i  = 0; i < S_i_.m; ++i){
  fprintf(debugFile,"\n%d %d\n", S_i_.row_number[i], S_i_.max_row_number[i]);
  for(unsigned int j = 0; j < S_i_.row_number[i]; ++j){
    fprintf(debugFile, "%d %lf %lf\n", S_i_.index[i][j], S_i_.value1[i][j], S_i_.value2[i][j]);
  }
}
fprintf(debugFile,"<<< SYSTEMMATRIX AC\n");
fflush(debugFile);
    WEM_pwl(&S_i_, waveletList, nodeList, elementTree, T_, nPatches, nLevels,
            SingleLayer, DoubleLayer, factor);
debugFile = fopen("debug.out", "a");
fprintf(debugFile,">>> SYSTEMMATRIX AW\n");
fprintf(debugFile, "%d %d\n", S_i_.m, S_i_.n);
for(unsigned int i  = 0; i < S_i_.m; ++i){
  fprintf(debugFile,"\n%d %d\n", S_i_.row_number[i], S_i_.max_row_number[i]);
  for(unsigned int j = 0; j < S_i_.row_number[i]; ++j){
    fprintf(debugFile, "%d %lf %lf\n", S_i_.index[i][j], S_i_.value1[i][j], S_i_.value2[i][j]);
  }
}
fprintf(debugFile,"<<< SYSTEMMATRIX AW\n");
fflush(debugFile);
    aposteriori1_ = postproc_pwl(&S_i_, waveletList, elementTree, nPatches, nLevels);
debugFile = fopen("debug.out", "a");
fprintf(debugFile,">>> SYSTEMMATRIX AP\n");
fprintf(debugFile, "%d %d\n", S_i_.m, S_i_.n);
for(unsigned int i  = 0; i < S_i_.m; ++i){
  fprintf(debugFile,"\n%d %d\n", S_i_.row_number[i], S_i_.max_row_number[i]);
  for(unsigned int j = 0; j < S_i_.row_number[i]; ++j){
    fprintf(debugFile, "%d %lf %lf\n", S_i_.index[i][j], S_i_.value1[i][j], S_i_.value2[i][j]);
  }
}
fprintf(debugFile,"<<< SYSTEMMATRIX AP\n");
fflush(debugFile);
}

void PWLSolver::constructSe()
{
    gf = greenOutside_; // sets the global pointer to pass GF to C code
    apriori2_ = compression_pwl(&S_e_, waveletList, elementTree, nPatches, nLevels,
                                nNodes);

FILE* debugFile = fopen("debug.out", "a");
fprintf(debugFile,">>> SYSTEMMATRIX AC\n");
fprintf(debugFile, "%d %d\n", S_i_.m, S_i_.n);
for(unsigned int i  = 0; i < S_i_.m; ++i){
  fprintf(debugFile,"\n%d %d\n", S_i_.row_number[i], S_i_.max_row_number[i]);
  for(unsigned int j = 0; j < S_i_.row_number[i]; ++j){
    fprintf(debugFile, "%d %lf %lf\n", S_i_.index[i][j], S_i_.value1[i][j], S_i_.value2[i][j]);
  }
}
fprintf(debugFile,"<<< SYSTEMMATRIX AC\n");
fflush(debugFile);
    WEM_pwl(&S_e_, waveletList, nodeList, elementTree, T_, nPatches, nLevels,
            SingleLayer, DoubleLayer, -2*M_PI);
debugFile = fopen("debug.out", "a");
fprintf(debugFile,">>> SYSTEMMATRIX AW\n");
fprintf(debugFile, "%d %d\n", S_i_.m, S_i_.n);
for(unsigned int i  = 0; i < S_i_.m; ++i){
  fprintf(debugFile,"\n%d %d\n", S_i_.row_number[i], S_i_.max_row_number[i]);
  for(unsigned int j = 0; j < S_i_.row_number[i]; ++j){
    fprintf(debugFile, "%d %lf %lf\n", S_i_.index[i][j], S_i_.value1[i][j], S_i_.value2[i][j]);
  }
}
fprintf(debugFile,"<<< SYSTEMMATRIX AW\n");
fflush(debugFile);
    aposteriori2_ = postproc_pwl(&S_e_, waveletList, elementTree, nPatches, nLevels);
debugFile = fopen("debug.out", "a");
fprintf(debugFile,">>> SYSTEMMATRIX AP\n");
fprintf(debugFile, "%d %d\n", S_i_.m, S_i_.n);
for(unsigned int i  = 0; i < S_i_.m; ++i){
  fprintf(debugFile,"\n%d %d\n", S_i_.row_number[i], S_i_.max_row_number[i]);
  for(unsigned int j = 0; j < S_i_.row_number[i]; ++j){
    fprintf(debugFile, "%d %lf %lf\n", S_i_.index[i][j], S_i_.value1[i][j], S_i_.value2[i][j]);
  }
}
fprintf(debugFile,"<<< SYSTEMMATRIX AP\n");
fflush(debugFile);
}

void PWLSolver::solveFirstKind(const Eigen::VectorXd & potential,
                               Eigen::VectorXd & charge)
{
    sparse G;
    double * rhs = 0;
    double * u = (double*) calloc(nNodes, sizeof(double));
    double * v = (double*) calloc(nNodes, sizeof(double));
    //next line is just a quick fix but i do not like it...
    Eigen::VectorXd pot = potential;
    WEMRHS2M_pwl(&rhs, waveletList, elementTree, T_, nPatches, nLevels, nNodes,
                 pot.data(), quadratureLevel_); // Transforms pot data to wavelet representation
    int iters = WEMPGMRES2_pwl(&S_i_, rhs, v, threshold, waveletList, elementList,
                               nPatches, nLevels); // v = A^{-1} * rhs
    init_sparse(&G, nNodes, nNodes, 10);
    single_scale_gram_pwl(&G, elementList, nPatches, nLevels);
    tdwtLin(v, elementList, nLevels, nPatches, nNodes);
    for (size_t i = 0; i < nNodes; i++) {
        for (size_t j = 0; j < G.row_number[i]; j++) {
            u[i] += G.value[i][j] * v[G.index[i][j]];
        }
    }
    dwtLin(u, elementList, nLevels, nPatches, nNodes);
    for (size_t i = 0; i < nNodes; i++) {
        rhs[i] += 4 * M_PI * u[i] / (epsilon -
                                     1); // assembling RHS equation (2.7) Computing, 2009, 86, 1-22
    }
    memset(u, 0, nNodes * sizeof(double));
    iters = WEMPCG_pwl(&S_i_, rhs, u, threshold, waveletList, elementList, nPatches,
                       nLevels);
    tdwtLin(u, elementList, nLevels, nPatches, nNodes);
    charge_pwl(u, charge.data(), elementList, T_, nPatches, nLevels);
    energy_pwl(u, pot.data(), elementList, T_, nPatches, nLevels);
    free(rhs);
    free(u);
    free(v);
    free_sparse(&G);
}

void PWLSolver::solveSecondKind(const Eigen::VectorXd & potential,
                                Eigen::VectorXd & charge)
{
    throw std::runtime_error("Second Kind (electric field) not yet implemented.");
}

void PWLSolver::solveFull(const Eigen::VectorXd & potential,
                          Eigen::VectorXd & charge)
{
    double * rhs = 0;
    double * u = (double*) calloc(nNodes, sizeof(double));
    double * v = (double*) calloc(nNodes, sizeof(double));
    //next line is just a quick fix but i do not like it...
    Eigen::VectorXd pot = potential;
    WEMRHS2M_pwl(&rhs, waveletList, elementTree, T_, nPatches, nLevels, nNodes,
                 pot.data(), quadratureLevel_); // Transforms pot data to wavelet representation
    int iters = WEMPCG_pwl(&S_i_, rhs, u, threshold, waveletList, elementList, nPatches,
                           nLevels);  /* u = V_i^(-1)*N_f */
    memset(rhs, 0, nNodes * sizeof(double));
    for (size_t i = 0; i < nNodes; i++) {
        /* rhs = V_e*u */
        for (size_t j = 0; j < S_e_.row_number[i]; j++) {
            rhs[i] += S_e_.value1[i][j] * u[S_e_.index[i][j]];
        }
    }
    iters = WEMPGMRES3_pwl(&S_i_, &S_e_, rhs, v, threshold, waveletList, elementList,
                           nPatches, nLevels); // v = A^{-1} * rhs
    for (size_t i = 0; i < nNodes; i++) {
        u[i] -= 4 * M_PI * v[i];      /* u = u - 4*pi*v */
    }
    tdwtLin(u, elementList, nLevels, nPatches, nNodes);
    charge_pwl(u, charge.data(), elementList, T_, nPatches, nLevels);
    energy_pwl(u, pot.data(), elementList, T_, nPatches, nLevels);
    free(rhs);
    free(u);
    free(v);
}

std::ostream & PWLSolver::printSolver(std::ostream & os)
{
    os << "Solver Type: Wavelet, piecewise linear functions";
    return os;
}
