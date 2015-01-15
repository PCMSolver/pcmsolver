/**
 * @file Mask.hpp
 *
 * @brief defines the fast wavelet transform masks for constant, linear wavelets
 */

#include "SparseMatrix.hpp"

// chooses according to m the correct mask
// constant - defines the mask T on level m
void mask_T11(SparseMatrix * T, unsigned int m);

// constant - defines the mask T on level m
void mask_T13(SparseMatrix* T, unsigned int m);

// linear - defines the mask T on level m, 2 vanishing moments
void mask_T22(SparseMatrix *T, unsigned int m, unsigned int M);

// linear - defines the mask T on level m, 4 vanishing moments
void mask_T24(SparseMatrix *T, unsigned int m, unsigned int M);

template <class AFClass> void dwtMask(SparseMatrix *T, SparseMatrix *L, unsigned int m, unsigned int M, AFClass *af);
