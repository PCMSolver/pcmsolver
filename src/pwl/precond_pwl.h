/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Weffc++"
#pragma clang diagnostic ignored "-Wextra"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdeprecated-register"
#pragma clang diagnostic ignored "-Wincompatible-pointer-types"
#pragma clang diagnostic ignored "-Wempty-body"
#endif

/* warning-disabler-end */

#ifndef PRECOND_PWL
#define PRECOND_PWL
/***************
 *  precond_pwl.h  *
 ***************/


/*================================================*
 *  Berechnet den Preconditioner und              *
 *  definiert dessen Anwendung auf einen Vektor.  *
 *================================================*/


void inv_A_times_x_pwl(sparse *A, double *x, unsigned int **F, unsigned int p, unsigned int M);
/* Loest das lineare Gleichungssystem T*A*T'*y = x,
   wobei der Nullvektor als Startvektor verwendet 
   und die Loesung y in x geschrieben wird. */


void single_scale_gram_pwl(sparse *G, unsigned int **F, unsigned int p, unsigned int M);
/* berechnet die Massenmatrix in der Einskalenbasis.
   Fuer die Quadratur wird die Mittelpunktsregel verwendet. */


void precond_pwl(double *a, double *b, sparse *G, wavelet_pwl * W, unsigned int **F, unsigned int p, unsigned int M);
/* Anwendung des Preconditioners auf den Vektor b, der
   NICHT veraendert wird. Das Ergebnis wird in a gespeichert */
#endif
/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

/* warning-disabler-end */

