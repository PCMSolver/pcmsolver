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

/****************
 *  Matlab_C.c  *
 ****************/


/*=================================================*
 *  Schnittstelle zwischen Matlab und C:           *
 *  Konvertiert die Patch- und Punkteliste         *
 *  vom Matlab-Format ins C-Format und umgekehrt.  *
 *=================================================*/


#include <stdlib.h>
#include "vector3.h"


void c2matlab_pointlist(P, PP, np)
/* Konvertiert die vector3-Punkteliste P in eine
   (3*np,1)-REAL-Matrix PP fuer Matlab und gibt
   den Speicherplatz von P frei. */
vector3 **P;
double *PP;
unsigned int np;
{
    unsigned int k;
    for (k = 0; k < np; k++) {
        PP[k] = (*P)[k].x;
        PP[k + np] = (*P)[k].y;
        PP[k + 2 * np] = (*P)[k].z;
    }
    free(*P);
    return;
}


vector3 *matlab2c_pointlist(PP, np)
/* Konvertiert die (3*np,1)-REAL-Punkteliste PP
   in eine vector3-Punkteliste fuer C */
double *PP;
unsigned int np;
{
    vector3 *P;
    unsigned int k;
    P = (vector3 *) calloc(np, sizeof(vector3));
    for (k = 0; k < np; k++) {
        P[k].x = PP[k];
        P[k].y = PP[k + np];
        P[k].z = PP[k + 2 * np];
    }
    return (P);
}


void c2matlab_patchlist(F, FF, nf)
/* Konvertiert die (nf,4)-(unsigned int)-Matrix F in eine
   (4*nf,1)-REAL-Matrix FF fuer Matlab, wobei die Indizes
   um eins erhoeht werden. Der Speicherplatz von F wird
   freigegeben. */
unsigned int ***F;
double *FF;
unsigned int nf;
{
    unsigned int k;
    for (k = 0; k < nf; k++) {
        FF[k] = (double) (*F)[k][0] + 1;
        FF[k + nf] = (double) (*F)[k][1] + 1;
        FF[k + 2 * nf] = (double) (*F)[k][2] + 1;
        FF[k + 3 * nf] = (double) (*F)[k][3] + 1;
        free((*F)[k]);
    }
    free(*F);
    return;
}


unsigned int **matlab2c_patchlist(FF, nf)
/* Konvertiert die (4*nf,1)-REAL-Patchliste aus Matlab
   in eine (nf,4)-(unsigned int)-Matrix F fuer C,
   wobei die Indizes um 1 reduziert werden. */
double *FF;
unsigned int nf;
{
    unsigned int **F;
    unsigned int k;
    F = (unsigned int **) calloc(nf, sizeof(unsigned int *));
    for (k = 0; k < nf; k++) {
        F[k] = (unsigned int *) calloc(4, sizeof(unsigned int));
        F[k][0] = (unsigned int) FF[k] - 1;
        F[k][1] = (unsigned int) FF[k + nf] - 1;
        F[k][2] = (unsigned int) FF[k + 2 * nf] - 1;
        F[k][3] = (unsigned int) FF[k + 3 * nf] - 1;
    }
    return (F);
}
/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

/* warning-disabler-end */

