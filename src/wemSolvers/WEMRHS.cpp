#include "GenericAnsatzFunction.hpp"
#include "GaussSquare.hpp"
#include "string.h"
#include "data.hpp"

#ifdef DEBUG2
#include <cstdio>
#endif

void WEMRHS1(double **rhs, GenericAnsatzFunction *af){
  unsigned int n = 1 << af->nLevels;
  double h = 1./n;
  Cubature *Q;
  double *c;
  double **y;
  Vector2 t;
  Vector3 n_t;
  double w;

  initGaussSquare(&Q,af->quadratureLevel_+1);
  y = (double**)malloc(af->totalSizeElementList*sizeof(double*)+af->noPhi*af->totalSizeElementList*sizeof(double));
  for(unsigned int i = 0; i < af->totalSizeElementList;++i) y[i] = (double*)(y+af->totalSizeElementList)+i*af->noPhi;
  c = (double*)malloc(af->noPhi*sizeof(double));
  (*rhs) = (double*) malloc(af->waveletList.sizeWaveletList*sizeof(double));

  // 1. Quadrature on fine level
  for(unsigned int i = af->noPatch*(n*n-1)/3; i <af->totalSizeElementList; ++i){
    memset(c,0,sizeof(double)*af->noPhi);
    for(unsigned int k = 0; k < Q[af->quadratureLevel_].noP; ++k){
      t.x = h*(af->elementTree.element[i].index_s+Q[af->quadratureLevel_].xi[k].x);
      t.y = h*(af->elementTree.element[i].index_t+Q[af->quadratureLevel_].xi[k].y);
      n_t = af->interCoeff->n_Chi(t, af->elementTree.element[i].patch);
      w = Q[af->quadratureLevel_].weight[k]*vector3Dot(df(af->interCoeff->Chi(t,af->elementTree.element[i].patch)),n_t);
      af->calculateCRHS(c,w,Q[af->quadratureLevel_].xi[k]);
    }
    for(unsigned int j = 0; j < af->noPhi; ++j){
      y[i][j] = h*c[j];
    }
  }
  //2. calculate coarser integrals from fine level
  for(int i = af->noPatch*(n*n-1)/3-1; i >=0; --i){
    af->calculateYRHS(y,i);
  }

  //3. set integrals (f,psi) together
  for(unsigned int i = 0; i < af->waveletList.sizeWaveletList;++i){
    w = 0;
    for(unsigned int j = 0; j < af->waveletList.W[i].noElements; ++j){
      for(unsigned int k = 0; k < af->noPhi;++k){
        w += y[af->waveletList.W[i].element[j]][k]*af->waveletList.W[i].weight[j*af->noPhi+k];
      }
    }
    (*rhs)[i] = w;
  }
  freeGaussSquare(&Q, af->quadratureLevel_+1);
  free(y);
  free(c);
  return;
}

void WEMRHS2(double **rhs, GenericAnsatzFunction *af){
  unsigned int n = 1 << af->nLevels;
  double h = 1./n;
  Cubature *Q;
  double *c;
  double **y;
  Vector2 t;
  double w;

  initGaussSquare(&Q,af->quadratureLevel_+1);
  y = (double**)malloc(af->totalSizeElementList*sizeof(double*)+af->noPhi*af->totalSizeElementList*sizeof(double));
  for(unsigned int i = 0; i < af->totalSizeElementList;++i) y[i] = (double*)(y+af->totalSizeElementList)+i*af->noPhi;
  c = (double*)malloc(af->noPhi*sizeof(double));
  (*rhs) = (double*) malloc(af->waveletList.sizeWaveletList*sizeof(double));

  // 1. Quadrature on fine level
  for(unsigned int i = af->noPatch*(n*n-1)/3; i <af->totalSizeElementList; ++i){
    memset(c,0,sizeof(double)*af->noPhi);
    for(unsigned int k = 0; k < Q[af->quadratureLevel_].noP; ++k){
      t.x = h*(af->elementTree.element[i].index_s+Q[af->quadratureLevel_].xi[k].x);
      t.y = h*(af->elementTree.element[i].index_t+Q[af->quadratureLevel_].xi[k].y);
      w = Q[af->quadratureLevel_].weight[k]*f(af->interCoeff->Chi(t,af->elementTree.element[i].patch));
      af->calculateCRHS(c,w,Q[af->quadratureLevel_].xi[k]);
    }
    for(unsigned int j = 0; j < af->noPhi; ++j){
      y[i][j] = h*c[j];
    }
  }
  //2. calculate coarser integrals from fine level
  for(int i = af->noPatch*(n*n-1)/3-1; i >=0; --i){
    af->calculateYRHS(y,i);
  }

  //3. set integrals (f,psi) together
  for(unsigned int i = 0; i < af->waveletList.sizeWaveletList; ++i){
    w = 0;
    for(unsigned int j = 0; j < af->waveletList.W[i].noElements; ++j){
      for(unsigned int k = 0; k < af->noPhi;++k){
        w += y[af->waveletList.W[i].element[j]][k]*af->waveletList.W[i].weight[j*af->noPhi+k];
      }
    }
    (*rhs)[i] = w;
  }
  freeGaussSquare(&Q, af->quadratureLevel_+1);
  free(y);
  free(c);
  return;
}

void WEMRHS2M(double **rhs, double *potential, GenericAnsatzFunction *af){
  unsigned int n = 1 << af->nLevels;
  double h = 1./n;
  Cubature *Q;
  double *c;
  double **y;
  Vector2 t;
  double w;
  unsigned int index;

  initGaussSquare(&Q,af->quadratureLevel_+1);
  y = (double**)malloc(af->totalSizeElementList*sizeof(double*)+af->noPhi*af->totalSizeElementList*sizeof(double));
  for(unsigned int i = 0; i < af->totalSizeElementList;++i) y[i] = (double*)(y+af->totalSizeElementList)+i*af->noPhi;
  c = (double*)malloc(af->noPhi*sizeof(double));
  (*rhs) = (double*) malloc(af->waveletList.sizeWaveletList*sizeof(double));

  // 1. Quadrature on fine level
  for(unsigned int i = af->noPatch*(n*n-1)/3; i <af->totalSizeElementList; ++i){
    memset(c,0,sizeof(double)*af->noPhi);
    for(unsigned int k = 0; k < Q[af->quadratureLevel_].noP; ++k){
      t.x = h*(af->elementTree.element[i].index_s+Q[af->quadratureLevel_].xi[k].x);
      t.y = h*(af->elementTree.element[i].index_t+Q[af->quadratureLevel_].xi[k].y);
      //w = Q[af->quadratureLevel_].weight[k]*f(af->interCoeff->Chi(t,af->elementTree.element[i].patch));
      index = (af->elementTree.element[i].patch*n*n + af->elementTree.element[i].index_t*n+af->elementTree.element[i].index_s);
      w = Q[af->quadratureLevel_].weight[k]*potential[index*Q[af->quadratureLevel].noP+k];
      af->calculateCRHS(c,w,Q[af->quadratureLevel_].xi[k]);
    }
    for(unsigned int j = 0; j < af->noPhi; ++j){
      y[i][j] = h*c[j];
    }
  }
  //2. calculate coarser integrals from fine level
  for(int i = af->noPatch*(n*n-1)/3-1; i >=0; --i){
    af->calculateYRHS(y,i);
  }

  //3. set integrals (f,psi) together
  for(unsigned int i = 0; i < af->waveletList.sizeWaveletList; ++i){
    w = 0;
    for(unsigned int j = 0; j < af->waveletList.W[i].noElements; ++j){
      for(unsigned int k = 0; k < af->noPhi;++k){
        w += y[af->waveletList.W[i].element[j]][k]*af->waveletList.W[i].weight[j*af->noPhi+k];
      }
    }
    (*rhs)[i] = w;
  }
  freeGaussSquare(&Q, af->quadratureLevel_+1);
  free(y);
  free(c);
  return;
}
