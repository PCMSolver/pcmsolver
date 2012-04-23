/**************
 *  Volume.c  *
 **************/


/*=======================================*
 * Berechnet das Volumen des Molekuels.  *
 *=======================================*/


#include <math.h>
#include <stdio.h>
#include "vector2.h"
#include "vector3.h"
#include "cubature.h"
#include "interpolate.h"
#include "gauss_square.h"
#include "volume.h"


double volume(F, T, p, m)
unsigned int **F;               /* Patchliste                                  */
vector3 ****T;                  /* Koeffizienten zur Oberflaecheninterpolation */
unsigned int p;                 /* Anzahl der Parametergebiete                 */
unsigned int m;                 /* Zahl der Level                              */
{
    unsigned int n = 1 << m;    /* n*n Patches pro Parametergebiet             */
    unsigned int i1, i2, i3;    /* Laufindizes fuer Ansatzfunktion             */
    cubature *Q;                /* Kubatur-Formeln                             */
    unsigned int g = 0;         /* Quadraturgrad                               */
    unsigned int l;             /* Laufindex fuer Quadratur                    */
    vector2 s;                  /* Linker, unterer Eckpunkt des Patches zi     */
    vector2 t;                  /* Auswertepunkte der Gaussquadratur           */
    double V;                   /* Energie                                     */
    double h = 1. / n;          /* Schrittweite                                */

/* Initialisierung */
    init_Gauss_Square(&Q, g + 1);       /* Kubatur-Formeln */

/* Berechnung des Fehlers */
    V = 0;
    for (i1 = 0; i1 < p; i1++) {
        s.y = 0;
        for (i2 = 0; i2 < n; i2++) {
            s.x = 0;
            for (i3 = 0; i3 < n; i3++) {        /* zeilenweise Durchnumerierung der Patches zi = (i1,i2,i3) */
                for (l = 0; l < Q[g].nop; l++) {
                    t = vector2_add(s, vector2_Smul(h, Q[g].xi[l]));
                    V += Q[g].w[l] * vector3_skalp(Chi(t, T[i1], m), n_Chi(t, T[i1], m));
                }
                s.x += h;
            }
            s.y += h;
        }
    }
    V *= h * h / 3;             /* correct scaling */

/* Datenausgabe */
    printf("Cavity's volume (> 0):           %g\n", V);
    free_Gauss_Square(&Q, g + 1);
    return (V);
}
