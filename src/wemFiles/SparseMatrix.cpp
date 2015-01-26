/**
 * @file SparseMatrix.cpp
 *
 * @brief sparse matrix implementation
 */
#include <string.h>
#include <stdlib.h>
#include "SparseMatrix.hpp"

/// copy a sparse matrix to another
void copySparse(SparseMatrix *from, SparseMatrix *to){
  unsigned int    i;

  to->m = from->m;
  to->n = from->n;
  if(from->value1){
    to->value1 = (double**) malloc(to->m*sizeof(double*));
    to->value2 = (double**) malloc(to->m*sizeof(double*));
    to->index = (unsigned int**) malloc(to->m*sizeof(unsigned int*));
    to->row_number = (unsigned int*) malloc(to->m*sizeof(unsigned int));
    to->max_row_number = (unsigned int*) malloc(to->m*sizeof(unsigned int));
    memcpy(to->row_number, from->row_number, to->m*sizeof(unsigned int));
    memcpy(to->max_row_number, from->max_row_number, to->m*sizeof(unsigned int));

    for (i=0; i<to->m; i++){
      to->index[i] = (unsigned int*) malloc((to->m+1)*sizeof(unsigned int));
      to->value1[i] = (double*) malloc((to->m)*sizeof(double));
      to->value2[i] = (double*) malloc((to->m)*sizeof(double));

      memcpy(to->value1[i], from->value1[i], (to->row_number[i])*sizeof(double));
      memcpy(to->value2[i], from->value2[i], to->row_number[i]*sizeof(double));
      memcpy(to->index[i], from->index[i], (to->row_number[i]+1)*sizeof(unsigned int));
    }
  }else{
    to->value1 = to->value2 = NULL;
    to->value = (double**) malloc(to->m*sizeof(double*));
    to->index = (unsigned int**) malloc(to->m*sizeof(unsigned int*));
    to->row_number = (unsigned int*) malloc(to->m*sizeof(unsigned int));
    to->max_row_number = (unsigned int*) malloc(to->m*sizeof(unsigned int));
    memcpy(to->row_number, from->row_number, to->m*sizeof(unsigned int));
    memcpy(to->max_row_number, from->max_row_number, to->m*sizeof(unsigned int));

    for (i=0; i<to->m; i++){
      to->index[i] = (unsigned int*) malloc((to->m+1)*sizeof(unsigned int));
      to->value[i] = (double*) malloc((to->m)*sizeof(double));

      memcpy(to->value[i], from->value[i], (to->row_number[i])*sizeof(double));
      memcpy(to->index[i], from->index[i], (to->row_number[i]+1)*sizeof(unsigned int));
    }
  }
}
/// binary search algorithm for element with index >= j
unsigned int searchSparse(unsigned int *array, unsigned int rn, unsigned int j) {
  unsigned int	mid, low, high;

  low  = 0;
  high = rn;
  while (low < high) {
    mid = (low+high)/2;
    if (array[mid] < j)
      low = mid+1;
    else if (array[mid] > j)
      high = mid;
    else return(mid);
  }
  // low == high!
  return(low);
}

/// memory move for use with insertion
void memmoveSparse(SparseMatrix *A, unsigned int i, unsigned int j) {
  unsigned int	rn, n_max;

  rn = A->row_number[i];

  // number of elements in row i is increased
  if (rn == A->max_row_number[i]){
    n_max = rn + 10;
    A->value[i] = (double*) realloc(A->value[i],n_max*sizeof(double));
    A->index[i] = (unsigned int*) realloc(A->index[i],(n_max+1)*sizeof(unsigned int));
    A->max_row_number[i] = n_max;
  }

  memmove(&A->value[i][j+1],&A->value[i][j],(rn-j)*sizeof(double));
  memmove(&A->index[i][j+1],&A->index[i][j],(rn+1-j)*sizeof(unsigned int));
  A->row_number[i]++;
  return;
}

/// memory move of index to use with insertion of pattern
void memmovePattern(SparseMatrix *A, unsigned int i, unsigned int j) {
  unsigned int	rn, n_max;

  rn = A->row_number[i];

  // number of elements in row i is increased
  if (rn == A->max_row_number[i]){
    n_max = rn + 10;
    A->index[i] = (unsigned int*) realloc(A->index[i],(n_max+1)*sizeof(unsigned int));
    A->max_row_number[i] = n_max;
  }

  memmove(&A->index[i][j+1],&A->index[i][j],(rn+1-j)*sizeof(unsigned int));
  A->row_number[i]++;
  return;
}

/// initialization routine for known size n,m with initial guess n_max non zero
/// elements per row
void initSparse(SparseMatrix *A, unsigned int m, unsigned int n, unsigned int n_max) {
  unsigned int    i;

  A->m = m;
  A->n = n;

  A->value = (double**) malloc(m*sizeof(double*));
  A->index = (unsigned int**) malloc(m*sizeof(unsigned int*));
  A->row_number = (unsigned int*) calloc(m,sizeof(unsigned int));
  A->max_row_number = (unsigned int*) malloc(m*sizeof(unsigned int));

  // insertion of dummy-element in index vector and assignment to n 
  for (i=0; i<m; i++){
    A->max_row_number[i] = n_max;
    A->value[i] = (double*) malloc(n_max*sizeof(double));
    A->index[i] = (unsigned int*) malloc((n_max+1)*sizeof(unsigned int));
    A->index[i][0] = n; // control element
  }
  return;
}

/// initialization of matrix pattern calculation - no value vector allocated
void initPattern(SparseMatrix *A, unsigned int m, unsigned int n, unsigned int n_max) {
  unsigned int    i;

  A->m = m;
  A->n = n;

  A->value = (double**) malloc(m*sizeof(double*));
  A->index = (unsigned int**) malloc(m*sizeof(unsigned int*));
  A->row_number = (unsigned int*) calloc(m,sizeof(unsigned int));
  A->max_row_number = (unsigned int*) malloc(m*sizeof(unsigned int));

  // insertion of dummy-element in index vector and assignment to n 
  for (i=0; i<m; i++){
    A->max_row_number[i] = n_max;
    A->index[i] = (unsigned int*) malloc((n_max+1)*sizeof(unsigned int));
    A->index[i][0] = n; // control element
  }

  return;
}

/// initialization of matrix pattern calculation for type 2 matrices - no value vector allocated
void initPattern2(SparseMatrix *A, unsigned int m, unsigned int n, unsigned int n_max) {
  unsigned int    i;

  A->m = m;
  A->n = n;

  A->value1 = (double**) malloc(m*sizeof(double*));
  A->value2 = (double**) malloc(m*sizeof(double*));
  A->index = (unsigned int**) malloc(m*sizeof(unsigned int*));
  A->row_number = (unsigned int*) calloc(m,sizeof(unsigned int));
  A->max_row_number = (unsigned int*) malloc(m*sizeof(unsigned int));

  // insertion of dummy-element in index vector and assignment to n 
  for (i=0; i<m; i++){
    A->max_row_number[i] = n_max;
    A->index[i] = (unsigned int*) malloc((n_max+1)*sizeof(unsigned int));
    A->index[i][0] = n; // control element
  }

  return;
}

/// add a pattern-element at row i and column j
void setPattern(SparseMatrix *A, unsigned int i, unsigned int j) {
  unsigned int	k;

  k = searchSparse(A->index[i],A->row_number[i],j);
  
  // entry not found
  if (A->index[i][k] != j){
    memmovePattern(A,i,k);
    A->index[i][k] = j;
  }

  return;
}

/// set value of row i column j
void setSparse(SparseMatrix *A,  unsigned int i, unsigned int j, double z){
  unsigned int	k;

  k = searchSparse(A->index[i],A->row_number[i],j);
  
  // entry not found
  if (A->index[i][k] != j){
    memmoveSparse(A,i,k);
    A->index[i][k] = j;
  }

  A->value[i][k] = z;
  return;
}

/// add in-place to row i, column j to value
void addSparse(SparseMatrix *A, unsigned int i, unsigned int j, double z) {
  unsigned int	k;

  if (z == 0) return;

  k = searchSparse(A->index[i],A->row_number[i],j);
  // element already in matrix
  if (A->index[i][k] == j) A->value[i][k] += z;
  // add element to matrix
  else {
    memmoveSparse(A,i,k);
    A->index[i][k] = j;
    A->value[i][k] = z;
  }
  return;
}

/// add in-place to row i, column j to value1 and value 2
void addSparse2(SparseMatrix *A, unsigned int i, unsigned int j, double z, double x) {
  unsigned int	k;

  if((z == 0) )
    if ((x == 0)) return;

  k = searchSparse(A->index[i],A->row_number[i],j);
  // element already in matrix
  if (A->index[i][k] == j){
    A->value1[i][k] += z;
    A->value2[i][k] += x;
  } else {
    // add element to matrix
    memmoveSparse(A,i,k);
    A->index[i][k] = j;
    A->value1[i][k] = z;
    A->value2[i][k] = x;
  }
  return;
}

/// retrieve value from row i column j
double getSparse(SparseMatrix *A, unsigned int i, unsigned int j) {
  unsigned int	k;

  k = searchSparse(A->index[i],A->row_number[i],j);
  // element found in matrix
  if (A->index[i][k] == j) return(A->value[i][k]);
  else return(0); // element not found
}

/// retrieve values from row i column j from value1 and value2
void getSparse2(SparseMatrix *A, unsigned int i, unsigned int j, double *z1, double *z2) {
  unsigned int	k;

  k = searchSparse(A->index[i],A->row_number[i],j);
  // element found in matrix
  if (A->index[i][k] == j){
    *z1 = (A->value1[i][k]);
    *z2 = (A->value2[i][k]);
  } else {
    // element not found
    *z1 = 0;
    *z2 = 0;
  }
}

/// count the non-zero elements of A
unsigned int nonZero(SparseMatrix *A){
  unsigned int    i, nz;

  nz = 0;
  for (i=0; i<A->m; i++) nz += A->row_number[i];

  return(nz);
}

/// reallocation of index and value vectors to proper size
void garbageCollect(SparseMatrix *A) {
  unsigned int	i, rn;

  for (i=0; i<A->m; i++){
    rn = A->row_number[i];
    A->value[i] = (double*) realloc(A->value[i],rn*sizeof(double));
    A->index[i] = (unsigned int*) realloc(A->index[i],(rn+1)*sizeof(unsigned int));
    A->max_row_number[i] = rn;
  }

  return;
}

/// allocation of value vector and reallocation of index to proper size after
/// calculation of non-zero element
void finishPattern(SparseMatrix *A) {
  unsigned int	i, rn;

  for (i=0; i<A->m; i++){
    rn = A->row_number[i];
    A->value[i] = (double*) calloc(rn,sizeof(double));
    A->index[i] = (unsigned int*) realloc(A->index[i],(rn+1)*sizeof(unsigned int));
    A->max_row_number[i] = rn;
    }

  return;
}

/// allocation of value vector1 and vector2 and reallocation of index to proper
/// size after calculation of non-zero element
void finishPattern2(SparseMatrix *A) {
  unsigned int	i, rn;

  for (i=0; i<A->m; i++){
    rn = A->row_number[i];
    A->value1[i] = (double*) calloc(rn,sizeof(double));
    A->value2[i] = (double*) calloc(rn,sizeof(double));
    A->index[i] = (unsigned int*) realloc(A->index[i],(rn+1)*sizeof(unsigned int));
    A->max_row_number[i] = rn;
  }

  return;
}

/// memory free function for type 1
void freeSparse(SparseMatrix *A) {
  unsigned int    i;

  for (i=0; i<A->m; i++){
    free(A->value[i]);
    free(A->index[i]);
  }

  free(A->value);
  free(A->index);
  free(A->row_number);
  free(A->max_row_number);

  return;
}

/// memory free function for type 2
void freeSparse2(SparseMatrix *A) {
  unsigned int    i;

  for (i=0; i<A->m; i++){
    free(A->value1[i]);
    free(A->value2[i]);
    free(A->index[i]);
  }

  free(A->value1);
  free(A->value2);
  free(A->index);
  free(A->row_number);
  free(A->max_row_number);

  return;
}
