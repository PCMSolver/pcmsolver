/**
 * @file Mask.hpp
 *
 * @brief defines the fast wavelet transform masks for constant, linear wavelets
 */

#include "SparseMatrix.hpp"

// chooses according to m the correct mask
template <class AFClass> void dwtMask(SparseMatrix *T, SparseMatrix *L, unsigned int m, unsigned int M, AFClass *af);
