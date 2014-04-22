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

/***********
 *  Pot.c  *
 ***********/


/*=============================================*
 *  Wertet das Potential in den Punkten R aus  *
 *  bzgl. des modifizierten Skalarproduktes.   *
 *=============================================*/


#include <math.h>
#include <stdlib.h>
#include "vector2.h"
#include "vector3.h"
#include "kern.h"
#include "cubature.h"
#include "interpolate.h"
#include "gauss_square.h"
#include "data.h"


#if !defined(pi)
#define pi 3.1415926535897932385
#endif


void pot(Pot, R, nr, u, T, p, m)
double **Pot;                   /* zu berechnendes Potential                   */
vector3 *R;                     /* Auswertepunkte des Potentials               */
unsigned int nr;                /* Anzahl der Auswertungspunkte                */
double *u;                      /* vorgegebene Dichte                          */
vector3 ****T;                  /* Koeffizienten zur Oberflaecheninterpolation */
unsigned int p;                 /* Anzahl der Patches                          */
unsigned int m;                 /* Zahl der Level                              */
{
    unsigned int n = 1 << m;    /* n*n Ansatzfunktionen pro Patch              */
    unsigned int i1, i2, i3;    /* Laufindizes fuer Ansatzfunktion             */
    unsigned int zi;            /* Zeilenindex hieraus: zi = i1*(n*n)+i2*n+i3  */
    unsigned int j;             /* Laufinidex fuer Auswertungspunkte           */
    cubature *Q;                /* Kubatur-Formeln                             */
    unsigned int g = 1;         /* Quadraturgrad                               */
    unsigned int k;             /* Laufindex fuer Quadratur                    */
    double w;                   /* Quadraturgewicht                            */
    vector2 s;                  /* Linker, unterer Eckpunkt des Patches zi     */
    vector2 t;                  /* Auswertepunkte der Gaussquadratur           */
    vector3 y;                  /* Auswertepunkt auf Gamma                     */
    vector3 n_y;                /* Normalableitung in obigem Punkt             */
    double h = 1. / n;          /* Schrittweite                                */

/* Initialisierung */
    init_Gauss_Square(&Q, g + 1);       /* Kubatur-Formeln */
    (*Pot) = (double *) calloc(nr, sizeof(double));

    zi = 0;
    for (i1 = 0; i1 < p; i1++) {
        s.y = 0;
        for (i2 = 0; i2 < n; i2++) {
            s.x = 0;
            for (i3 = 0; i3 < n; i3++) {        /* zeilenweise Durchnumerierung der Patches zi = (i1,i2,i3) */
                for (k = 0; k < Q[g].nop; k++) {
                    t = vector2_add(s, vector2_Smul(h, Q[g].xi[k]));
                    y = Chi(t, T[i1], m);
                    n_y = n_Chi(t, T[i1], m);
                    w = u[zi] / h;
                    for (j = 0; j < nr; j++) {  /* interior: Single - Double */
                        /* exterior: Double - Single */
                        (*Pot)[j] += Q[g].w[k] * (-w * SingleLayer(R[j], y) + DoubleLayer(R[j], y, n_y) * f(y));
                    }
                }
                s.x += h;
                zi++;
            }
            s.y += h;
        }
    }

/* Zum Schluss normieren */
    for (j = 0; j < nr; j++)
        (*Pot)[j] *= (h * h) / (4 * pi);

/* und den Speicherplatz wieder freigeben */
    free_Gauss_Square(&Q, g + 1);
    return;
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

