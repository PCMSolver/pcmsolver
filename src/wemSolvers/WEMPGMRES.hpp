/**
 * @file WEMPGMRES.hpp
 *
 * @brief implementation of several GMRES solvers
 */
#ifndef WEMPGMRES_HPP
#define WEMPGMRES_HPP
class GenericAnsatzFunction;

/**
 * @brief{ GMRES solver for linear system A2'*x = b
 * preconditioning through diagonal scaling}
 *
 * @param A the sparse matrix
 * @param b the RHS
 * @param x start value, will be replaced by solution
 * @param epsi precision
 * @param af the AnsatzFunction class for acces to wavelets, elements and
 * other constants
 */
unsigned int WEMPGMRES1(SparseMatrix *A, double *b, double *x, double epsi, 
	GenericAnsatzFunction *af);

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
unsigned int WEMPGMRES2(SparseMatrix *A, double *b, double *x, double epsi, 
	GenericAnsatzFunction *af);

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
unsigned int WEMPGMRES3(SparseMatrix *A, SparseMatrix *B, double *rhs, double *x, double epsi, 
	GenericAnsatzFunction *af);
#endif
