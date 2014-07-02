/**
 * @file Mask.hpp
 *
 * @brief defines the fast wavelet transform masks for constant with 3
 * vanishin moments and linear wavelets with 4 vanishing moments and
 * double points at patch edges
 */

//#include "Mask.hpp"
//#include "LinAnsatzFunction.hpp"
//#include "ConAnsatzFunction.hpp"
class LinAnsatzFunction;
class ConAnsatzFunction;

#include "SparseMatrix.hpp"
unsigned int     minLevel;

// constant - defines the mask T on level m
void mask_T11(SparseMatrix * T, unsigned int m);

// constant - defines the mask T on level m
void mask_T13(SparseMatrix* T, unsigned int m);

// linear - defines the mask T on level m, 2 vanishing moments
void mask_T22(SparseMatrix *T, unsigned int m, unsigned int M);

// linear - defines the mask T on level m, 4 vanishing moments
void mask_T24(SparseMatrix *T, unsigned int m, unsigned int M);

// choose the correct mask according to level
template <class D> void dwtMask(SparseMatrix *T, SparseMatrix *L, unsigned int m, unsigned int M, D *aF);

