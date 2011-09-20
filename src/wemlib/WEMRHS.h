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
