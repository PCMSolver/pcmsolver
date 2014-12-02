#include "GenericAnsatzFunction.hpp"
#include "GaussSquare.hpp"
#include "string.h"
//#include "data.hpp" - to replace with information from GreenFunctions
#include <cstdio>

#ifdef DEBUG2
#include <cstdio>
#endif

void WEMRHS2M(double **rhs, double *potential, GenericAnsatzFunction *af){
  unsigned int n = 1 << af->nLevels;
  double h = 1./n;
  Cubature *Q;
  double *c;
  double **y;
  Vector2 t;
  double w;
  unsigned int index;
#ifdef DEBUG
  FILE* debugFile = fopen("debug.out", "a");
  fprintf(debugFile,">>> POTENTIAL %d\n", af->quadratureLevel_);
#endif
  initGaussSquare(&Q,af->quadratureLevel_+1);
  y = (double**)malloc(af->elementTree.totalSizeElementList*sizeof(double*)+af->noPhi*af->elementTree.totalSizeElementList*sizeof(double));
  for(unsigned int i = 0; i < af->elementTree.totalSizeElementList;++i) y[i] = (double*)(y+af->elementTree.totalSizeElementList)+i*af->noPhi;
  c = (double*)malloc(af->noPhi*sizeof(double));
  (*rhs) = (double*) malloc(af->waveletList.sizeWaveletList*sizeof(double));

  // 1. Quadrature on fine level
  for(unsigned int i = af->nPatches*(n*n-1)/3; i <af->elementTree.totalSizeElementList; ++i){
    memset(c,0,sizeof(double)*af->noPhi);
    for(unsigned int k = 0; k < Q[af->quadratureLevel_].noP; ++k){
      t.x = h*(af->elementTree.element[i].index_s+Q[af->quadratureLevel_].xi[k].x);
      t.y = h*(af->elementTree.element[i].index_t+Q[af->quadratureLevel_].xi[k].y);
      //w = Q[af->quadratureLevel_].weight[k]*f(af->interCoeff->Chi(t,af->elementTree.element[i].patch));
      index = (af->elementTree.element[i].patch*n*n + af->elementTree.element[i].index_t*n+af->elementTree.element[i].index_s);
      w = Q[af->quadratureLevel_].weight[k]*potential[index*Q[af->quadratureLevel_].noP+k];
#ifdef DEBUG
      fprintf(debugFile,"%d %lf\n",index*Q[af->quadratureLevel_].noP+k, potential[index*Q[af->quadratureLevel_].noP+k]);
#endif
      af->calculateCRHS(c,w,Q[af->quadratureLevel_].xi[k]);
    }
    for(unsigned int j = 0; j < af->noPhi; ++j){
      y[i][j] = h*c[j];
    }
  }
#ifdef DEBUG
  fprintf(debugFile,"<<< POTENTIAL\n");
  fclose(debugFile);
#endif
  //2. calculate coarser integrals from fine level
  for(int i = af->nPatches*(n*n-1)/3-1; i >=0; --i){
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

const unsigned int nAtoms = 12;

/// position of charges: vectors (x,y,z)
const Vector3 x[12] = {
  { 5.274/0.52917721092,   1.999/0.52917721092,  -8.568/0.52917721092 },
  { 6.627/0.52917721092,   2.018/0.52917721092,  -8.209/0.52917721092 },
  { 7.366/0.52917721092,   0.829/0.52917721092,  -8.202/0.52917721092 },
  { 6.752/0.52917721092,  -0.379/0.52917721092,  -8.554/0.52917721092 },
  { 5.399/0.52917721092,  -0.398/0.52917721092,  -8.912/0.52917721092 },
  { 4.660/0.52917721092,   0.791/0.52917721092,  -8.919/0.52917721092 },
  { 4.704/0.52917721092,   2.916/0.52917721092,  -8.573/0.52917721092 },
  { 7.101/0.52917721092,   2.950/0.52917721092,  -7.938/0.52917721092 },
  { 8.410/0.52917721092,   0.844/0.52917721092,  -7.926/0.52917721092 },
  { 7.322/0.52917721092,  -1.296/0.52917721092,  -8.548/0.52917721092 },
  { 4.925/0.52917721092,  -1.330/0.52917721092,  -9.183/0.52917721092 },
  { 3.616/0.52917721092,   0.776/0.52917721092,  -9.196/0.52917721092 }
};

/// charges: double values
const double alpha[12]={6,6,6, 6,6,6, 1,1,1, 1,1,1}; 

double f(Vector3 a){
    double		c = 0;
    unsigned int	i;
    for (i=0; i<nAtoms; i++) c += alpha[i]/vector3Norm(vector3Sub(a, x[i]));
    return(c);
}

Vector3 df(Vector3 a){
    unsigned int	i;
    Vector3		c, r, v;
    c.x = c.y = c.z = 0;
    for (i=0; i<nAtoms; i++) {
        r = vector3Sub(a, x[i]);
        v = vector3SMul(alpha[i]/pow(vector3Norm(r),3),r);
        c.x -= v.x;
        c.y -= v.y;
        c.z -= v.z;
    }
    return(c);
}

void WEMRHS2M_test(double **rhs, double *potential, GenericAnsatzFunction *af){
  unsigned int n = 1 << af->nLevels;
  double h = 1./n;
  Cubature *Q;
  double *c;
  double **y;
  Vector2 t;
  double w;
  unsigned int index;
#ifdef DEBUG
  FILE* debugFile = fopen("debug.out", "a");
  fprintf(debugFile,">>> POTENTIAL %d\n", af->quadratureLevel_);
#endif
  initGaussSquare(&Q,af->quadratureLevel_+1);
  y = (double**)malloc(af->elementTree.totalSizeElementList*sizeof(double*)+af->noPhi*af->elementTree.totalSizeElementList*sizeof(double));
  for(unsigned int i = 0; i < af->elementTree.totalSizeElementList;++i) y[i] = (double*)(y+af->elementTree.totalSizeElementList)+i*af->noPhi;
  c = (double*)malloc(af->noPhi*sizeof(double));
  (*rhs) = (double*) malloc(af->waveletList.sizeWaveletList*sizeof(double));

  // 1. Quadrature on fine level
  for(unsigned int i = af->nPatches*(n*n-1)/3; i <af->elementTree.totalSizeElementList; ++i){
    memset(c,0,sizeof(double)*af->noPhi);
    for(unsigned int k = 0; k < Q[af->quadratureLevel_].noP; ++k){
      t.x = h*(af->elementTree.element[i].index_s+Q[af->quadratureLevel_].xi[k].x);
      t.y = h*(af->elementTree.element[i].index_t+Q[af->quadratureLevel_].xi[k].y);
      //w = Q[af->quadratureLevel_].weight[k]*f(af->interCoeff->Chi(t,af->elementTree.element[i].patch));
      index = (af->elementTree.element[i].patch*n*n + af->elementTree.element[i].index_t*n+af->elementTree.element[i].index_s);
      w = Q[af->quadratureLevel_].weight[k]*vector3Dot(df(af->interCoeff->Chi(t,af->elementTree.element[i].patch)),af->interCoeff->n_Chi(t, af->elementTree.element[i].patch));
#ifdef DEBUG
      fprintf(debugFile,"%d %lf\n",index*Q[af->quadratureLevel_].noP+k, potential[index*Q[af->quadratureLevel_].noP+k]);
#endif
      af->calculateCRHS(c,w,Q[af->quadratureLevel_].xi[k]);
    }
    for(unsigned int j = 0; j < af->noPhi; ++j){
      y[i][j] = h*c[j];
    }
  }
#ifdef DEBUG
  fprintf(debugFile,"<<< POTENTIAL\n");
  fclose(debugFile);
#endif
  //2. calculate coarser integrals from fine level
  for(int i = af->nPatches*(n*n-1)/3-1; i >=0; --i){
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
