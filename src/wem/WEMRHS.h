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

#ifndef WEMRHS
#define WEMRHS
/**************
 *  WEMRHS.h  *
 **************/


/*===================================*
 *  Bestimmt die rechte Seite im     *
 *  Wavelet-Galerkin-Verfahren fuer  *
 *  stueckweise konstante Wavelets.  *
 *===================================*/


void WEMRHS1(double **rhs, wavelet *W, element *E, vector3 ****T, unsigned int p, unsigned int M);
/* testet die Neumann-Daten des gegebenen Potentials */


void WEMRHS2(double **rhs, wavelet *W, element *E, vector3 ****T, unsigned int p, unsigned int M);

void WEMRHS2M(double **rhs, wavelet *W, element *E, vector3 ****T, unsigned int p, unsigned int M, double *potential, unsigned int g);
/* testet die Dirichlet-Daten des gegebenen Potentials */
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

