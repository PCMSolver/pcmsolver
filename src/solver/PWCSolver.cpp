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
#include "PWCSolver.hpp"

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
#include "energy.h"
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

void PWCSolver::initPointers()
{
    elementTree = NULL;
    waveletList = NULL;
}


/*PWCSolver::PWCSolver(const Section & solver) : WEMSolver(solver)
  {
	initPointers();
	setSolverType("Wavelet");
}*/

PWCSolver::~PWCSolver()
{
    if(elementTree != NULL) free_elementlist(&elementTree, nPatches,nLevels);
    if(waveletList != NULL) free_waveletlist(&waveletList, nPatches,nLevels);
    if(T_ != NULL)          free_interpolate(&T_,nPatches,nLevels);
}

void PWCSolver::initInterpolation()
{
    init_interpolate(&T_, pointList, nPatches, nLevels);
    nNodes = gennet(&nodeList, &elementList, pointList, nPatches, nLevels);
#ifdef DEBUG
    FILE* debugFile = fopen("debug.out","w");
    fclose(debugFile);
    debugFile = fopen("debug.out","a");	
    fprintf(debugFile,">>> PPOINTLIST\n");
    for(unsigned int m = 0; m<nNodes; ++m){
      fprintf(debugFile,"%lf %lf %lf\n",(nodeList)[m].x, (nodeList)[m].y, (nodeList)[m].z);
    }
    fprintf(debugFile,"<<< PPOINTLIST\n");
    fprintf(debugFile,">>> PELEMENTLIST\n");
    unsigned int nf = nPatches*(1<<nLevels)*(1<<nLevels);
    for(unsigned int m = 0; m<nf; ++m){
      fprintf(debugFile,"%d %d %d %d\n",(elementList)[m][0], (elementList)[m][1], (elementList)[m][2], (elementList)[m][3]);
    }
    fprintf(debugFile,"<<< PELEMENTLIST\n");
    fclose(debugFile);
#endif
}

void PWCSolver::constructWavelets()
{
    generate_elementlist(&elementTree, nodeList, elementList, nPatches, nLevels);
    generate_waveletlist(&waveletList, elementTree, nPatches, nLevels);
    set_quadrature_level(waveletList, elementTree, nPatches, nLevels);
    simplify_waveletlist(waveletList, elementTree, nPatches, nLevels);
    complete_elementlist(waveletList, elementTree, nPatches, nLevels);
#ifdef DEBUG
    unsigned int	N = 1 << nLevels;	
    unsigned int	ne;
    ne = nPatches*(4*N*N-1)/3;
    FILE* debugFile = fopen("debug.out","a");
    fprintf(debugFile,">>> WAVELET_TREE_SIMPLIFY\n");
    for(unsigned int m = 0; m<nNodes; ++m){
            fprintf(debugFile,"%d %d %d %d\n", m, waveletList[m].level, waveletList[m].element_number, 4);
            for(unsigned int i1 = 0; i1< waveletList[m].element_number;++i1){
                    fprintf(debugFile,"%d %lf ", waveletList[m].element[i1], waveletList[m].weight[i1]);
            }
            for(unsigned int i1 = 0; i1< 4;++i1){
                    fprintf(debugFile,"%d ", waveletList[m].son[i1]);
            }
            fprintf(debugFile,"\n");
    }
    fprintf(debugFile,"<<< WAVELET_TREE_SIMPLIFY\n");
    fflush(debugFile);
    fprintf(debugFile,">>> HIERARCHICAL_ELEMENT_TREE\n");
    for(unsigned int m = 0; m<ne; ++m){
        fprintf(debugFile,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf\n", elementTree[m].patch, elementTree[m].level, elementTree[m].index_s, elementTree[m].index_t, elementTree[m].midpoint.x, elementTree[m].midpoint.y, elementTree[m].midpoint.z, elementTree[m].radius, nodeList[elementTree[m].vertex[0]].x, nodeList[elementTree[m].vertex[0]].y, nodeList[elementTree[m].vertex[0]].z);
        for(unsigned int i1 = 0; i1< elementTree[m].wavelet_number;++i1)
            fprintf(debugFile,"%d ", elementTree[m].wavelet[i1]);
        fprintf(debugFile,"\n");
    }
    fprintf(debugFile,"<<< HIERARCHICAL_ELEMENT_TREE\n");
    fclose(debugFile);
#endif
}

void PWCSolver::constructSi()
{
    double factor = 0;
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
    apriori1_ = compression(&S_i_, waveletList, elementTree, nPatches, nLevels);
#ifdef DEBUG
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
#endif
    WEM(&S_i_, waveletList, elementTree, T_, nPatches, nLevels, SingleLayer, DoubleLayer,
        factor);
#ifdef DEBUG
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
#endif

    aposteriori1_ = postproc(&S_i_, waveletList, elementTree, nPatches, nLevels);
#ifdef DEBUG
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
#endif

}

void PWCSolver::constructSe()
{
    gf = greenOutside_; // sets the global pointer to pass GF to C code
    apriori2_ = compression(&S_e_, waveletList, elementTree, nPatches, nLevels);
    WEM(&S_e_, waveletList, elementTree, T_, nPatches, nLevels, SingleLayer, DoubleLayer,
        -2*M_PI);
    aposteriori2_ = postproc(&S_e_, waveletList, elementTree, nPatches, nLevels);
}

void PWCSolver::solveFirstKind(const Eigen::VectorXd & potential,
                               Eigen::VectorXd & charge)
{
    double *rhs;
    double *u = (double*) calloc(nFunctions, sizeof(double));
    double * pot = const_cast<double *>(potential.data());
    double * chg = charge.data();
    double epsilon = greenOutside_->epsilon();
    WEMRHS2M(&rhs, waveletList, elementTree, T_, nPatches, nLevels, pot,
             quadratureLevel_);
#ifdef DEBUG
    FILE *debugFile = fopen("debug.out", "a");
    fprintf(debugFile,">>> WEMRHS1\n");
    for(unsigned int i = 0; i < nNodes; ++i){
      fprintf(debugFile,"%d %lf\n",i, rhs[i]);
    }
    fprintf(debugFile,"<<< WEMRHS1\n");
    fflush(debugFile);
#endif

    int iter = WEMPGMRES2(&S_i_, rhs, u, threshold, nPatches, nLevels);
#ifdef DEBUG
    fprintf(debugFile,">>> WEMPGMRES1\n");
    for(unsigned int i = 0; i < nNodes; ++i){
      fprintf(debugFile,"%d %lf\n",i, u[i]);
    }
    fprintf(debugFile,"<<< WEMPGMRES1\n");
    fflush(debugFile);
#endif

    tdwtKon(u, nLevels, nFunctions);
    dwtKon(u, nLevels, nFunctions);
    for (size_t i = 0; i < nFunctions; ++i) {
        rhs[i] += 4 * M_PI * u[i] / (epsilon - 1);
    }
#ifdef DEBUG
    fprintf(debugFile,">>> WEMRHS2\n");
    for(unsigned int i = 0; i < nNodes; ++i){
      fprintf(debugFile,"%d %lf\n",i, rhs[i]);
    }
    fprintf(debugFile,"<<< WEMRHS2\n");
    fflush(debugFile);
#endif

    memset(u, 0, nFunctions * sizeof(double));
    iter = WEMPCG(&S_i_, rhs, u, threshold, nPatches, nLevels);
#ifdef DEBUG
    fprintf(debugFile,">>> WEMPCG %g\n",threshold);
    for(unsigned int i = 0; i <  nNodes; ++i){
      fprintf(debugFile,"%d %.10lf\n",i, u[i]);
    }
    fprintf(debugFile,"<<< WEMPCG\n");
    fflush(debugFile);
#endif

    tdwtKon(u, nLevels, nFunctions);
    energy_ext(u, pot, elementList, T_, nPatches, nLevels);
 #ifdef DEBUG
    unsigned int	ne;		/* Anzahl der Elemente                        */
    unsigned int	N = 1 << nLevels;	/* N*N Elemente pro Patch auf dem Level M     */
    ne = nPatches*(4*N*N-1)/3;			/* Anzahl der Elemente */
    fprintf(debugFile,">>> CHARGE %g\n",threshold);
    for(unsigned int i = 0; i <  ne; ++i){
      fprintf(debugFile,"%d %.10lf\n",i, charge[i]);
    }
    fprintf(debugFile,"<<< CHARGE\n");
    fflush(debugFile);
    fclose(debugFile);
#endif
   charge_ext(u, chg, elementList, T_, nPatches, nLevels);
    free(rhs);
    free(u);
}

void PWCSolver::solveSecondKind(const Eigen::VectorXd & potential,
                                Eigen::VectorXd & charge)
{
    throw std::runtime_error("Second Kind not yet implemented"); //careful to the double layer sign when implementing it....
}

void PWCSolver::solveFull(const Eigen::VectorXd & potential,
                          Eigen::VectorXd & charge)
{
    double *rhs;
    double *u = (double*) calloc(nFunctions, sizeof(double));
    double *v = (double*) calloc(nFunctions, sizeof(double));
    //next line is just a quick fix to avoid problems with const but i do not like it...
    double * pot = const_cast<double *>(potential.data());
    WEMRHS2M(&rhs, waveletList, elementTree, T_, nPatches, nLevels, pot,
             quadratureLevel_);
    int iters = WEMPCG(&S_i_, rhs, u, threshold, nPatches, nLevels);
    memset(rhs, 0, nFunctions*sizeof(double));
    for(unsigned int i = 0; i < nFunctions; i++) {
        for(unsigned int j = 0; j < S_e_.row_number[i]; j++) {
            rhs[i] += S_e_.value1[i][j] * u[S_e_.index[i][j]];
        }
    }
    iters = WEMPGMRES3(&S_i_, &S_e_, rhs, v, threshold,
                       nPatches, nLevels);
    for(unsigned int i = 0; i < nFunctions; i++) {
        u[i] -= 4*M_PI*v[i];
    }
    tdwtKon(u, nLevels, nFunctions);
    energy_ext(u, pot, elementList, T_, nPatches, nLevels);
    charge_ext(u, charge.data(), elementList, T_, nPatches, nLevels);
    free(rhs);
    free(u);
    free(v);
}

std::ostream & PWCSolver::printSolver(std::ostream & os)
{
    os << "Solver Type: Wavelet, piecewise constant functions";
    return os;
}
