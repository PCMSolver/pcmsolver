/**
 * @file SparseMatrix.hpp
 *
 * @brief sparse matrix implementation
 */

#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP

typedef struct {
  unsigned int    m, n; ///< the size of the matrix
  unsigned int    *row_number, *max_row_number; ///< pointer to the index and value vectors where the specific row starts and the maximum row number
  unsigned int	   **index; ///< column vector
  double    **value;  ///< value pointer used for simple matrices
  double    **value1; ///< 1st value pointer used for double matrices
  double    **value2; ///< 2nd value pointer used for double matrices
} SparseMatrix;

/// initialization routine - each row n_max elements allocated
void initSparse(SparseMatrix *A, 
	unsigned int m, unsigned int n, unsigned int n_max);

/// initialization of matrix pattern
void initPattern(SparseMatrix *A, 
	unsigned int m, unsigned int n, unsigned int n_max);
/// initialization of matrix pattern - each row n_max elements allocated
void initPattern2(SparseMatrix *A, 
	unsigned int m, unsigned int n, unsigned int n_max);
    
/// binary search algorithm for element with index >= j
unsigned int searchSparse(unsigned int *array, unsigned int rn, unsigned int j);

/// add a pattern-element at row i and column j
void setPattern(SparseMatrix *A, unsigned int i, unsigned int j);

/// set value of row i column j
void setSparse(SparseMatrix *A, unsigned int i, unsigned int j, double z);

/// add in-place to row i, column j to value
void addSparse(SparseMatrix *A, unsigned int i, unsigned int j, double z);
/// add in-place to row i, column j to value1 and value 2
void addSparse2(SparseMatrix *A, unsigned int i, unsigned int j, double z, double x);

/// retrieve value from row i column j
double getSparse(SparseMatrix *A, unsigned int i, unsigned int j);
/// retrieve values from row i column j from value1 and value2
void getSparse2(SparseMatrix *A, unsigned int i, unsigned int j, double* v1, double *v2);

/// count the non-zero elements of A
unsigned int nonZero(SparseMatrix *A);

/// reallocation of index and value vectors to proper size
void garbageCollect(SparseMatrix *A);

/// allocation of value vector and reallocation of index to proper size after
/// calculation of non-zero element
void finishPattern(SparseMatrix *A);
/// allocation of value vector1 and vector2 and reallocation of index to proper
/// size after calculation of non-zero element
void finishPattern2(SparseMatrix *A);

/// memory free function for type 1
void freeSparse(SparseMatrix *A);
/// memory free function for type 2
void freeSparse2(SparseMatrix *A);
/// copy a sparse matrix to another
void copySparse(SparseMatrix *A1, SparseMatrix *A2);
#endif
