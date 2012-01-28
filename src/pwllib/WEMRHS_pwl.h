#ifndef WEMRHS_pwl
#define WEMRHS_pwl
/**************
 *  WEMRHS_pwl.h  *
 **************/


/*===================================*
 *  Bestimmt die rechte Seite im     *
 *  Wavelet-Galerkin-Verfahren fuer  *
 *  stueckweise lineare Wavelets.    *
 *===================================*/


void WEMRHS_pwl1(double **rhs, wavelet *W, element *E, vector3 ****T, unsigned int p, unsigned int M, unsigned int np);
/* testet die Neumann-Daten des gegebenen Potentials */

void WEMRHS_pwl2(double **rhs, wavelet *W, element *E, vector3 ****T, unsigned int p, unsigned int M, unsigned int np);
/* testet die Dirichlet-Daten des gegebenen Potentials */

void WEMRHS_pwl2M(double **rhs, wavelet *W, element *E, vector3 ****T, unsigned int p, unsigned int M, unsigned int np, double *potential, unsigned int g);
#endif
