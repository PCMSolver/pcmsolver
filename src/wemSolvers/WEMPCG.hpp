/*
 * @file WEMPCG.hpp
 *
 * @brief the CG solver
 */
#ifndef CG_HPP
#define CG_HPP
unsigned int WEMPCG(SparseMatrix* A, double* b, double* x, double epsi, GenericAnsatzFunction* af);
#endif
