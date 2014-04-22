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

/*********************
 *  init_exterior.c  *
 *********************/


/*========================================================*
 *  Berechnet die Auswertepunkte auf der Einheitssphaere  *
 *========================================================*/


#include <math.h>
#include <stdlib.h>
#include "vector3.h"
#include "init_points.h"


#if !defined pi
#define pi 3.1415926535897932385
#endif


double R = 1.5;                 /* Radius */


unsigned int init_points(Q, m)
vector3 **Q;
unsigned int m;
{
    double *a;
    unsigned int n = 1 << m;
    unsigned int nq;
    unsigned int i, j;
    unsigned int zq;
    double h = 1. / n;

/* Initialisierung */
    nq = 6 * n * n;             /* Anzahl der Auswertungspunkte */
    (*Q) = (vector3 *) malloc(nq * sizeof(vector3));
    a = (double *) malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
        a[i] = h * (i + 0.5);

/* Berechne Auswertepunkte */
    zq = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            (*Q)[zq++] = vector3_make(+1, 2 * a[i] - 1, 2 * a[j] - 1);
            (*Q)[zq++] = vector3_make(-1, 2 * a[i] - 1, 2 * a[j] - 1);
            (*Q)[zq++] = vector3_make(2 * a[i] - 1, +1, 2 * a[j] - 1);
            (*Q)[zq++] = vector3_make(2 * a[i] - 1, -1, 2 * a[j] - 1);
            (*Q)[zq++] = vector3_make(2 * a[i] - 1, 2 * a[j] - 1, +1);
            (*Q)[zq++] = vector3_make(2 * a[i] - 1, 2 * a[j] - 1, -1);
        }
    }
    for (i = 0; i < nq; i++)
        (*Q)[i] = vector3_Smul(R / vector3_norm((*Q)[i]), (*Q)[i]);

/* Speicherplatz wieder freigeben */
    free(a);
    return (nq);
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

