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
 *  Topology.c  *
 ****************/


/*====================================*
 *  Erstellt die topologischen Daten  *
 *====================================*/


#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "vector2.h"
#include "vector3.h"
#include "topology_pwl.h"


void init_grid_pwl(vector3 **P, unsigned int ***F, vector3 ***U, unsigned int p, unsigned int m, unsigned int *np, unsigned int *nf);


void refine_grid_pwl(vector3 **P, unsigned int ***F, vector3 ***U, unsigned int p, unsigned int m, unsigned int M, unsigned int *np, unsigned int *nf);


void init_grid_pwl(P, F, U, p, m, np, nf)
/* Erstellt die Punkte- und Indexliste fuer Level 0. */
vector3 **P;                    /* Zeiger auf die Punkteliste          */
unsigned int ***F;              /* Zeiger auf die lokale Basisliste    */
vector3 ***U;                   /* Gitterpunkte                        */
unsigned int p;                 /* Anzahl der Parametergebiete         */
unsigned int m;                 /* 2^m*2^m Patches pro Parametergebiet */
unsigned int *np, *nf;          /* sizeof(P) bzw. sizeof(F)            */
{
    unsigned int n = 1 << m;    /* n*n Patches pro Parametergebiet     */
    unsigned int i;             /* Laufindizes fuer Parametergebiet    */

/* Speicherplatz allokieren */
    (*P) = (vector3 *) malloc(4 * p * sizeof(vector3));
    (*F) = (unsigned int **) malloc(p * sizeof(unsigned int *));

/* Eckpunkte der Parametergebiete bestimmen */
    *np = 4 * p;
    *nf = p;
    for (i = 0; i < p; i++) {
        /* Bestimme die 4 Eckpunkte des Patches */
        (*P)[4 * i] = U[i][0][0];
        (*P)[4 * i + 1] = U[i][0][n];
        (*P)[4 * i + 2] = U[i][n][0];
        (*P)[4 * i + 3] = U[i][n][n];

        (*F)[i] = (unsigned int *) malloc(4 * sizeof(unsigned int));
        (*F)[i][0] = 4 * i;
        (*F)[i][1] = 4 * i + 1;
        (*F)[i][2] = 4 * i + 3;
        (*F)[i][3] = 4 * i + 2;
    }
    return;
}


void refine_grid_pwl(P, F, U, p, m, M, np, nf)
/* Erstellt die Punkte- und Indexliste fuer alle zusaetzlichen
   Gitterpunkte des Level m. */
vector3 **P;                    /* Zeiger auf die Punkteliste    */
unsigned int ***F;              /* Zeiger auf die Indexliste     */
vector3 ***U;                   /* Gitterpunkte                  */
unsigned int p;                 /* Anzahl der Parametergebiete   */
unsigned int m;                 /* aktuelles Level               */
unsigned int M;                 /* hoechstes Level               */
unsigned int *np, *nf;          /* sizeof(P) bzw. sizeof(F)      */
{
    unsigned int n = 1 << (m - 1);      /* n = 2^(m-1)                   */
    unsigned int S = 1 << (M - m);      /* Schrittweite zum nexten Punkt */
    signed int i1, i2, i3;      /* Laufindizes                   */
    unsigned int fz;            /* Patchzaehler                  */

/* Speicherplatz allokieren */
    (*P) = (vector3 *) realloc(*P, 4 * (*np) * sizeof(vector3));
    (*F) = (unsigned int **) realloc(*F, 4 * (*nf) * sizeof(unsigned int *));
    for (i1 = *nf; i1 < 4 * (*nf); i1++)
        (*F)[i1] = (unsigned int *) calloc(4, sizeof(unsigned int));

/* Kopiere altes Gitter in neues Gitter */
    fz = (p - 1) * n * n;
    for (i1 = p - 1; i1 >= 0; i1--) {
        for (i2 = n - 1; i2 >= 0; i2--) {
            for (i3 = n - 1; i3 >= 0; i3--) {
                (*F)[4 * fz + 2 * n * (2 * i2 + 1) + 2 * i3][3] = (*F)[fz + n * i2 + i3][3];
                (*F)[4 * fz + 2 * n * (2 * i2 + 1) + 2 * i3 + 1][2] = (*F)[fz + n * i2 + i3][2];
                (*F)[4 * fz + 4 * n * i2 + 2 * i3 + 1][1] = (*F)[fz + n * i2 + i3][1];
                (*F)[4 * fz + 4 * n * i2 + 2 * i3][0] = (*F)[fz + n * i2 + i3][0];
            }
        }
        fz -= n * n;
    }

/* Bestimme Verfeinerung */
    fz = 0;
    *nf *= 4;
    for (i1 = 0; i1 < p; i1++) {
        for (i2 = 0; i2 <= 2 * n; i2++) {
            for (i3 = 0; i3 <= 2 * n; i3++) {
                if ((i2 % 2 == 1) || (i3 % 2 == 1)) {
                    (*P)[*np] = U[i1][i2 * S][i3 * S];
                    if ((i2 < 2 * n) && (i3 < 2 * n))
                        (*F)[fz + 2 * n * i2 + i3][0] = *np;
                    if ((i2 < 2 * n) && (0 < i3))
                        (*F)[fz + 2 * n * i2 + i3 - 1][1] = *np;
                    if ((0 < i2) && (0 < i3))
                        (*F)[fz + 2 * n * (i2 - 1) + i3 - 1][2] = *np;
                    if ((0 < i2) && (i3 < 2 * n))
                        (*F)[fz + 2 * n * (i2 - 1) + i3][3] = *np;
                    (*np)++;
                }
            }
        }
        fz += 4 * n * n;
    }

/* Speicherplatz wieder freigeben */
    return;
}


unsigned int gennet_pwl(P, F, U, p, M)
/* Erstellt die Punkte- und Patchliste in hierarchischer Weise
   und liefert als Funktionsergebnis die Laenge der Punkteliste,
   die vom Geschlecht der Oberflaeche abhaengig ist. */
vector3 **P;                    /* Zeiger auf die Punkteliste          */
unsigned int ***F;              /* Zeiger auf die Patchliste           */
vector3 ***U;                   /* Gitterpunkte                        */
unsigned int p;                 /* Anzahl der Patches                  */
unsigned int M;                 /* 2^M*2^M Patches pro Parametergebiet */
{
    unsigned int m;             /* Laufindex fuer das Level            */
    unsigned int np;            /* Laenge von P                        */
    unsigned int nf;            /* Laenge von F                        */

    init_grid_pwl(P, F, U, p, M, &np, &nf);
    for (m = 1; m <= M; m++)
        refine_grid_pwl(P, F, U, p, m, M, &np, &nf);
    return (np);
}


void free_patchlist_pwl(F, nf)
/* gibt den Speicherplatz der (nf,4)-(unsigned int)-Patchliste F frei */
unsigned int ***F;
unsigned int nf;
{
    unsigned int k;
    for (k = 0; k < nf; k++)
        free((*F)[k]);
    free(*F);
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

