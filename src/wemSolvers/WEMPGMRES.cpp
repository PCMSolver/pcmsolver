/**
 * @file WEMPGMRES.cpp
 *
 * @brief implementation of several GMRES solvers
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "SparseMatrix.hpp"
#include "Vector3.hpp"
#include "WEMPGMRES.hpp"
#include "GenericAnsatzFunction.hpp"
#ifdef DEBUG2
#include <cstdio>
#endif
const unsigned int 	maxiter = 100; // restart parameter

/*=================*
 *  Hauptprogramm  *
 *=================*/
 
/**
 * @brief GMRES solver for linear system A2'*x = b
 * preconditioning through diagonal scaling
 *
 * @param A the sparse matrix
 * @param b the RHS
 * @param x start value, will be replaced by solution
 * @param epsi precision
 * @param af the AnsatzFunction class for acces to wavelets, elements and
 * other constants
 */
unsigned int WEMPGMRES1(SparseMatrix *A, double *b, double *x, double epsi, GenericAnsatzFunction * /* af */) {
  unsigned int	k, l;
  double  	rn, t;
  double  	**H, **V, *D, *v, *f, *u, *c, *s, *y;

  V = (double**) malloc((maxiter+1)*sizeof(double*));
  H = (double**) malloc((maxiter+1)*sizeof(double*));

  for (unsigned int i=0; i<=maxiter; ++i) {
    V[i] = (double*) malloc(A->n*sizeof(double));
    H[i] = (double*) malloc(maxiter*sizeof(double));
  }
   
  v = (double*) malloc(A->n*sizeof(double));
  u = (double*) malloc((maxiter+1)*sizeof(double));
  f = (double*) malloc((maxiter+1)*sizeof(double));
  c = (double*) malloc(maxiter*sizeof(double));
  s = (double*) malloc(maxiter*sizeof(double));
  y = (double*) malloc(maxiter*sizeof(double));

  // preconditioning: D = diag(A)^(-1)
  D = (double*) malloc(A->n*sizeof(double));
  for (unsigned int i=0; i<A->n; ++i) {
    getSparse2(A,i,i,&rn,&t);
    D[i] = 1/fabs(t);
  }
  rn = 1;
  k = 0;

  // outer iteration
  for (l=0; rn>epsi; l++){  
    // initial residual V(0,:) = D*(b - A*x), rn = sqrt(V(0,:)*V(0,:)') 
    rn = 0;
    memcpy(V[0],b,A->n*sizeof(double));
    for (unsigned int i=0; i<A->n; ++i){
      for (unsigned int j=0; j<A->row_number[i]; ++j){
        V[0][A->index[i][j]] -= A->value2[i][j] * x[i];
      }
    }
    for (unsigned int i=0; i<A->n; ++i){
      V[0][i] *= D[i];
      rn += V[0][i]*V[0][i];
    }
    rn = sqrt(rn);

    // V(0,:) /= rn
    t = 1/rn;
    for (unsigned int i=0; i<A->n; ++i) V[0][i] *= t;
    u[0] = rn;
   
    // inner iteration
    for (k=0; (k<maxiter) && (rn>epsi); k++){  
      // matrix vector product v = D*(A*V(k,:))
      memset(v,0,A->n*sizeof(double));
      for (unsigned int i=0; i<A->n; ++i){
        for (unsigned int j=0; j<A->row_number[i]; ++j){
          v[A->index[i][j]] += A->value2[i][j] * V[k][i];
        }
      }
      for (unsigned int i=0; i<A->n; ++i) v[i] *= D[i];	 

      // calculate scalps f = V(0:k,:)*v
      memset(f,0,(k+1)*sizeof(double));
      for (unsigned int i=0; i<=k; ++i){
        for (unsigned int j=0; j<A->n; ++j) f[i] += V[i][j] * v[j];
      }

      // (k+1)-th basisvector v = v - V(0:k,:)*f and rn = sqrt(v'*v)
      rn = 0;
      for (unsigned int i=0; i<A->n; ++i){
        for (unsigned int j=0; j<=k; ++j) v[i] -= V[j][i] * f[j];
        rn += v[i]*v[i];
      }
      rn = sqrt(rn);

      // V(k+1,:) = v' / rn
      t = 1/rn;
      for (unsigned int i=0; i<A->n; ++i) V[k+1][i] = v[i] * t;

      // calculate Givens rotation
      t = f[0];
      for (unsigned int i=0; i<k; ++i){
        f[i] = c[i] * t + s[i] * f[i+1];
        t = c[i] * f[i+1] - s[i] * t;
      }
      f[k] = sqrt(t*t + rn*rn);
      c[k] = t  / f[k];
      s[k] = rn / f[k];

      // H(0:k,k) = f
      for (unsigned int i=0; i<=k; ++i) H[i][k] = f[i];
      u[k+1] = - s[k] * u[k];
      u[k]   = c[k] * u[k];
      rn = fabs(u[k+1]);
    }// end inner iteration
   
    // solving triangular linear system y = H(0:k-1,0:k-1) \ u(0:k-1) 
    for (int i=k-1; i>=0; --i){
      t = 0;
      for (int j=k-1; j>i; --j) t += H[i][j] * y[j];
      y[i] = (u[i] - t) / H[i][i];
    }
   
    // new initial guess x = x + V(0:k-1,:)'*y 
    for (unsigned int i=0; i<A->n; ++i){
      for (unsigned int j=0; j<k; ++j) x[i] += V[j][i] * y[j];
    }
  }
  // free memory
  for (unsigned int i=0; i<=maxiter; ++i) {
    free(V[i]);
    free(H[i]);
  }
  free(V);
  free(H);
  free(D);
  free(v);
  free(f);
  free(u); 
  free(c);
  free(s);
  free(y);

  return((l-1)*maxiter+k);
}

/**
 * @brief{ GMRES solver for linear system A2*x = b
 * preconditioning through diagonal scaling}
 *
 * @param A the sparse matrix
 * @param b the RHS
 * @param x start value, will be replaced by solution
 * @param epsi precision
 * @param af the AnsatzFunction class for acces to wavelets, elements and
 * other constants
 */
unsigned int WEMPGMRES2(SparseMatrix *A, double *b, double *x, double epsi, GenericAnsatzFunction * /* af */) {
  unsigned int	k, l;
  double  	rn, t;
  double  	**H, **V, *D, *v, *f, *u, *c, *s, *y;

  V = (double**) malloc((maxiter+1)*sizeof(double*));
  H = (double**) malloc((maxiter+1)*sizeof(double*));

  for (unsigned int i=0; i<=maxiter; ++i) {
    V[i] = (double*) malloc(A->n*sizeof(double));
    H[i] = (double*) malloc(maxiter*sizeof(double));
  }
   
  v = (double*) malloc(A->n*sizeof(double));
  u = (double*) malloc((maxiter+1)*sizeof(double));
  f = (double*) malloc((maxiter+1)*sizeof(double));
  c = (double*) malloc(maxiter*sizeof(double));
  s = (double*) malloc(maxiter*sizeof(double));
  y = (double*) malloc(maxiter*sizeof(double));

  // Preconditioning: D = diag(A)^(-1)
  D = (double*) malloc(A->n*sizeof(double));
  for (unsigned int i=0; i<A->n; ++i){
    getSparse2(A,i,i,&rn,&t);
    D[i] = 1/fabs(t);
  }
  rn = 1;
  k = 0;

  // outer iteration
  for (l=0; rn>epsi; l++){  
    // initial residual: V(0,:) = D*(b - A*x), rn = sqrt(V(0,:)*V(0,:)')
    rn = 0;
    memcpy(V[0],b,A->n*sizeof(double));
    for (unsigned int i=0; i<A->n; ++i){
      for (unsigned int j=0; j<A->row_number[i]; ++j){
        V[0][i] -= A->value2[i][j] * x[A->index[i][j]];
      }
      V[0][i] *= D[i];
      rn += V[0][i]*V[0][i];
    }
    rn = sqrt(rn);

    // V(0,:) /= rn
    t = 1/rn;
    for (unsigned int i=0; i<A->n; ++i) V[0][i] *= t;
    u[0] = rn;

    // inner iteration
    for (k=0; (k<maxiter) && (rn>epsi); ++k){  
      // matrix vector product v = D*(A*V(k,:))
      memset(v,0,A->n*sizeof(double));
      for (unsigned int i=0; i<A->n; ++i){
        for (unsigned int j=0; j<A->row_number[i]; ++j){
          v[i] += A->value2[i][j] * V[k][A->index[i][j]];
        }
        v[i] *= D[i];
      }

      // calculate scalps? f = V(0:k,:)*v
      memset(f,0,(k+1)*sizeof(double));
      for (unsigned int i=0; i<=k; ++i){
        for (unsigned int j=0; j<A->n; ++j) f[i] += V[i][j] * v[j];
      }

      // (k+1)-th basisvector v = v - V(0:k,:)*f and rn = sqrt(v'*D*v)
      rn = 0;
      for (unsigned int i=0; i<A->n; ++i){
        for (unsigned int j=0; j<=k; ++j) v[i] -= V[j][i] * f[j];
        rn += v[i]*v[i];
      }
      rn = sqrt(rn);

      // V(k+1,:) = v' / rn
      t = 1/rn;
      for (unsigned int i=0; i<A->n; ++i) V[k+1][i] = v[i] * t;

      // calculate Givens rotation
      t = f[0];
      for (unsigned int i=0; i<k; ++i){
        f[i] = c[i] * t + s[i] * f[i+1];
        t = c[i] * f[i+1] - s[i] * t;
      }

      f[k] = sqrt(t*t + rn*rn);
      c[k] = t  / f[k];
      s[k] = rn / f[k];

      // H(0:k,k) = f
      for (unsigned int i=0; i<=k; ++i) H[i][k] = f[i];

      u[k+1] = - s[k] * u[k];
      u[k]   = c[k] * u[k];
      rn = fabs(u[k+1]);
    }// end inner iteration

    // solve triangular linear system y = H(0:k-1,0:k-1) \ u(0:k-1)
    for (int i=k-1; i>=0; --i){
      t = 0;
      for (int j=k-1; j>i; --j) t += H[i][j] * y[j];
      y[i] = (u[i] - t) / H[i][i];
    }

    // new initial guess x = x + V(0:k-1,:)'*y
    for (unsigned int i=0; i<A->n; ++i){
      for (unsigned int j=0; j<k; ++j) x[i] += V[j][i] * y[j];
    }
  }

  // free memory
  for (unsigned int i=0; i<=maxiter; ++i){
    free(V[i]);
    free(H[i]);
  }

  free(V);
  free(H);
  free(D);
  free(v);
  free(f);
  free(u); 
  free(c);
  free(s);
  free(y);

  return((l-1)*maxiter+k);;
}

/**
 * @brief{ GMRES solver for linear system (B1*G^(-1)*A2'-B2*G^(-1)*A1)*x = rhs
 * preconditioning through wavelet scaling}
 *
 * @param A the sparse matrix
 * @param B the sparse matrix
 * @param rhs the RHS
 * @param x start value, will be replaced by solution
 * @param epsi precision
 * @param af the AnsatzFunction class for acces to wavelets, elements and
 * other constants
 */
unsigned int WEMPGMRES3(SparseMatrix *A, SparseMatrix *B, double *rhs, double *x, double epsi, GenericAnsatzFunction *af){
  unsigned int	k, l;
  double  	rn, t;
  double  	**H, **V, *v, *f, *u, *c, *s, *y;

  // calculate gramian
  af->createGram(A->n,10);
  
  // allocate memory
  V = (double**) malloc((maxiter+1)*sizeof(double*));
  H = (double**) malloc((maxiter+1)*sizeof(double*));

  for (unsigned int i=0; i<=maxiter; ++i){
    V[i] = (double*) malloc(A->n*sizeof(double));
    H[i] = (double*) malloc(maxiter*sizeof(double));
  }
   
  v = (double*) malloc(A->n*sizeof(double));
  u = (double*) malloc((maxiter+1)*sizeof(double));
  f = (double*) malloc((maxiter+1)*sizeof(double));
  c = (double*) malloc(maxiter*sizeof(double));
  s = (double*) malloc(maxiter*sizeof(double));
  y = (double*) malloc(maxiter*sizeof(double));

  rn = 1;
  k = 0;

  // outer iteration 
  for (l=0; rn>epsi; l++){  
    // inital residual V(0,:) = rhs - (B1*G^(-1)*A2'-B2*G^(-1)*A1)*x
    memset(v,0,A->n*sizeof(double));
    memcpy(V[0],rhs,A->n*sizeof(double));
    for (unsigned int i=0; i<A->n; ++i){
      for (unsigned int j=0; j<A->row_number[i]; ++j){
        v[A->index[i][j]] += A->value2[i][j] * x[i];
      }
    }
    af->inv_A_times_x(v);
    for (unsigned int i=0; i<B->n; ++i){
      for (unsigned int j=0; j<B->row_number[i]; ++j){
        V[0][i] -= B->value1[i][j] * v[B->index[i][j]];
      }
    }
    memset(v,0,A->n*sizeof(double));
    for (unsigned int i=0; i<A->n; ++i){
      for (unsigned int j=0; j<A->row_number[i]; ++j){
        v[i] += A->value1[i][j] * x[A->index[i][j]];
      }
    }
    af->inv_A_times_x(v);
    for (unsigned int i=0; i<B->n; ++i){
      for (unsigned int j=0; j<B->row_number[i]; ++j){
        V[0][i] += B->value2[i][j] * v[B->index[i][j]];
      }
    }	 
   	 
    // v = precond(V(0,:)) and rn = sqrt(v'*v) 	 
    rn = 0;
    af->precond(v,V[0]);
    for (unsigned int i=0; i<A->n; ++i) rn += v[i]*v[i];
    rn = sqrt(rn);

    // v = V(0,:)/rn
    t = 1./rn;
    for (unsigned int i=0; i<A->n; ++i) V[0][i] = t * v[i];
    u[0] = rn;

    // inner iteration
    for (k=0; (k<maxiter) && (rn>epsi); k++){  
      // matrix vector product V[k+1][i] = (B1*G^(-1)*A2'-B2*G^(-1)*A1)*V(k,:)
      memset(v,0,A->n*sizeof(double));
      memset(V[k+1],0,A->n*sizeof(double));
      for (unsigned int i=0; i<A->n; ++i){
        for (unsigned int j=0; j<A->row_number[i]; ++j){
          v[A->index[i][j]] += A->value2[i][j] * V[k][i];
        }
      }
      af->inv_A_times_x(v);
      for (unsigned int i=0; i<B->n; ++i){
        for (unsigned int j=0; j<B->row_number[i]; ++j){
          V[k+1][i] += B->value1[i][j] * v[B->index[i][j]];
        }
      }
      memset(v,0,A->n*sizeof(double));
      for (unsigned int i=0; i<A->n; ++i){
        for (unsigned int j=0; j<A->row_number[i]; ++j){
          v[i] += A->value1[i][j] * V[k][A->index[i][j]];
        }
      }
      af->inv_A_times_x(v);
      for (unsigned int i=0; i<B->n; ++i){
        for (unsigned int j=0; j<B->row_number[i]; ++j){
          V[k+1][i] -= B->value2[i][j] * v[B->index[i][j]];
        }
      }

      // v = precond(V(k+1,:)
      af->precond(v,V[k+1]);
      
      // calculate scalps? f = V(0:k,:)*v
      memset(f,0,(k+1)*sizeof(double));
      for (unsigned int i=0; i<=k; ++i){
        for (unsigned int j=0; j<A->n; ++j) f[i] += V[i][j] * v[j];
      }

      // (k+1)-th basisvector v = v - V(0:k,:)*f and rn = sqrt(v'*v)
      rn = 0;
      for (unsigned int i=0; i<A->n; ++i){
        for (unsigned int j=0; j<=k; ++j) v[i] -= V[j][i] * f[j];
        rn += v[i]*v[i];
      }
      rn = sqrt(rn);

      // V(k+1,:) = v' / rn
      t = 1./rn;
      for (unsigned int i=0; i<A->n; ++i) V[k+1][i] = v[i] * t;

      // calculate Givens rotation
      t = f[0];
      for (unsigned int i=0; i<k; ++i){
        f[i] = c[i] * t + s[i] * f[i+1];
        t = c[i] * f[i+1] - s[i] * t;
      }
      f[k] = sqrt(t*t + rn*rn);
      c[k] = t  / f[k];
      s[k] = rn / f[k];

      // H(0:k,k) = f 
      for (unsigned int i=0; i<=k; ++i) H[i][k] = f[i];

      u[k+1] = - s[k] * u[k];
      u[k]   = c[k] * u[k];
      rn = fabs(u[k+1]);
    } // end inner iteration

    // solve triangular linear system y = H(0:k-1,0:k-1) \ u(0:k-1)
    for (int i=k-1; i>=0; --i){
      t = 0;
      for (int j=k-1; j>i; --j) t += H[i][j] * y[j];
      y[i] = (u[i] - t) / H[i][i];
    }

    // new initial guess x = x + V(0:k-1,:)'*y
    for (unsigned int i=0; i<A->n; ++i){
      for (unsigned int j=0; j<k; ++j) x[i] += V[j][i] * y[j];
    }
  }

  // free memory
  af->freeGram();
  for (unsigned int i=0; i<=maxiter; ++i){
    free(V[i]);
    free(H[i]);
  }
  free(V);
  free(H);
  free(v);
  free(f);
  free(u); 
  free(c);
  free(s);
  free(y);
  return((l-1)*maxiter+k);
}
