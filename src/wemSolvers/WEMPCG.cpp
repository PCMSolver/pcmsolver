#include "SparseMatrix.hpp"
#include "GenericAnsatzFunction.hpp"
#include <cstdlib>
#include <string.h>

/*
 * @file WEMPCG.cpp
 *
 * @brief implementation of the preconditioned CG solver
 */

unsigned int WEMPCG(SparseMatrix* A, double* b, double* x, double epsi, GenericAnsatzFunction* af){
  double *r, *d, *z, *Ad;
  double u,v,omg;
  unsigned int k;

  // allocate memory
  r = (double*) malloc(A->n*sizeof(double));
  d = (double*) calloc(A->n,sizeof(double));
  z = (double*) malloc(A->n*sizeof(double));
  Ad = (double*) malloc(A->n*sizeof(double));

  // r = b-A*x
  memcpy(r,b,A->n*sizeof(double));
  for(unsigned int i = 0; i < A->n; ++i){
    for(unsigned int j = 0; j < A->row_number[i]; ++j){
      r[i] -= A->value1[i][j]*x[A->index[i][j]];
    }
  }

  // eventually create gram Matrix
  af->createGram(A->n, 10);

  // d = precond(r,M)
  af->precond(d,r);

  // u = (r,d)
  u = 0;
  for(unsigned int i = 0; i < A->n; ++i) u +=r[i]*d[i];

  // Iteration
  for(k = 0; sqrt(u) > epsi; k++){
    // Ad = A*d
    memset(Ad,0,A->n*sizeof(double));
    for(unsigned int i = 0; i < A->n; ++i){
      for(unsigned int j = 0; j < A->row_number[i]; ++j){
        Ad[i] += A->value1[i][j]*d[A->index[i][j]];
      }
    }

    // omg = u/(d,Ad) v = u
    omg = 0;
    for(unsigned int i = 0; i < A->n; ++i) omg += d[i]*Ad[i];
    omg = u/omg;
    v = u;

    // x = x+omg
    for(unsigned int i = 0; i < A->n; ++i){
      x[i] += omg *d[i];
      r[i] -= omg*Ad[i];
    }

    // z = precond(r,M)
    af->precond(z,r);

    // u = (r,z)
    u = 0;
    for(unsigned int i = 0; i < A->n; ++i) u+= r[i]*z[i];

    // d = z+u/v*d
    omg = u/v;
    for(unsigned int i = 0; i < A->n; ++i) d[i] = z[i]+omg*d[i];
  }

  af->freeGram();
  free(Ad);
  free(r);
  free(d);
  free(z);
  return k;
}
