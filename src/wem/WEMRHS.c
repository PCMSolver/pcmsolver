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

/**************
 *  WEMRHS.c  *
 **************/


/*===========================================================================*
 *  Bestimmt die rechte Seite im Wavelet-Galerkin-Verfahren fuer             *
 *  stueckweise konstante Wavelets bzgl. des modifizierten Skalarproduktes.  *
 *===========================================================================*/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "constants.h"
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "cubature.h"
#include "gauss_square.h"
#include "interpolate.h"
#include "WEMRHS.h"
//#include "data.h"
#include "molecule.h"

void WEMRHS1(rhs, W, E, T, p, M)
/* testet die Neumann-Daten des gegebenen Potentials */
double **rhs;                   /* zu berechnende rechte Seite                */
wavelet *W;                     /* Waveletliste                               */
element *E;                     /* hierarchische Elementliste                 */
vector3 ****T;                  /* Oberflaecheninterpolation                  */
unsigned int p;                 /* Anzahl der Patches                         */
unsigned int M;                 /* Zahl der Level                             */
{
    unsigned int N = 1 << M;    /* N*N Elemente pro Patch auf dem Level M     */
    unsigned int nw;            /* Anzahl der Wavelets                        */
    unsigned int ne;            /* Anzahl der Elemente                        */
    signed int i, j;            /* Laufindizes durch die Wavelet/Elementliste */
    double h;                   /* Schrittweite                               */
    cubature *Q;                /* Kubatur-Formeln                            */
    unsigned int g = 1;         /* Quadratur-Grad auf dem feinsten Level      */
    double c;                   /* Wert der Integrals                         */
    vector2 t;                  /* Auswertepunkt der Gauss-Quadratur          */
    vector3 n_t;                /* Auswertung der Normale an den obigen Punkt */
    unsigned int k;             /* Laufindex fuer die Quadratur               */
    double *y;                  /* Array mit den Integralen (f,phi_m)         */

    nw = p * N * N;             /* Anzahl der Wavelets */
    ne = p * (4 * N * N - 1) / 3;       /* Anzahl der Elemente */
    init_Gauss_Square(&Q, g + 1);       /* Kubatur-Formeln     */
    y = (double *) malloc(ne * sizeof(double));
    (*rhs) = (double *) malloc(nw * sizeof(double));

/* 1. Quadratur auf dem feinsten Level */
    h = 1. / N;
    for (i = ne - nw; i < ne; i++) {
        c = 0;
        for (k = 0; k < Q[g].nop; k++) {
            t.x = h * (E[i].index_s + Q[g].xi[k].x);
            t.y = h * (E[i].index_t + Q[g].xi[k].y);
            n_t = n_Chi(t, T[E[i].patch], M);
            c += Q[g].w[k] * vector3_skalp(field(Chi(t, T[E[i].patch], M)), n_t);
        }
        y[i] = h * c;
    }

/* 2. berechne Integrale der goeberen Level aus denen der feineren */
    for (i = ne - nw - 1; i >= 0; i--) {
        y[i] = 0.5 * (y[E[i].son[0]] + y[E[i].son[1]] + y[E[i].son[2]] + y[E[i].son[3]]);
    }

/* 3. setze Integrale (f,psi) zusammen */
    for (i = 0; i < nw; i++) {
        c = 0;
        for (j = 0; j < W[i].element_number; j++)
            c += y[W[i].element[j]] * W[i].weight[j];
        (*rhs)[i] = c;
    }

/* Speicherplatz wieder freigeben */
    free_Gauss_Square(&Q, g + 1);
    free(y);
    return;
}


void WEMRHS2(rhs, W, E, T, p, M)
/* testet die Dirichlet-Daten des gegebenen Potentials */
double **rhs;                   /* zu berechnende rechte Seite                */
wavelet *W;                     /* Waveletliste                               */
element *E;                     /* hierarchische Elementliste                 */
vector3 ****T;                  /* Oberflaecheninterpolation                  */
unsigned int p;                 /* Anzahl der Patches                         */
unsigned int M;                 /* Zahl der Level                             */
{
    unsigned int N = 1 << M;    /* N*N Elemente pro Patch auf dem Level M     */
    unsigned int nw;            /* Anzahl der Wavelets                        */
    unsigned int ne;            /* Anzahl der Elemente                        */
    signed int i, j;            /* Laufindizes durch die Wavelet/Elementliste */
    double h;                   /* Schrittweite                               */
    cubature *Q;                /* Kubatur-Formeln                            */
    unsigned int g = 1;         /* Quadratur-Grad auf dem feinsten Level      */
    double c;                   /* Wert der Integrals                         */
    vector2 t;                  /* Auswertepunkt der Gauss-Quadratur          */
    unsigned int k;             /* Laufindex fuer die Quadratur               */
    double *y;                  /* Array mit den Integralen (f,phi_m)         */

    nw = p * N * N;             /* Anzahl der Wavelets */
    ne = p * (4 * N * N - 1) / 3;       /* Anzahl der Elemente */
    init_Gauss_Square(&Q, g + 1);       /* Kubatur-Formeln     */
    y = (double *) malloc(ne * sizeof(double));
    (*rhs) = (double *) malloc(nw * sizeof(double));

/* 1. Quadratur auf dem feinsten Level */
    h = 1. / N;
    for (i = ne - nw; i < ne; i++) {
        c = 0;
        for (k = 0; k < Q[g].nop; k++) {
            t.x = h * (E[i].index_s + Q[g].xi[k].x);
            t.y = h * (E[i].index_t + Q[g].xi[k].y);
            c += Q[g].w[k] * potmol(Chi(t, T[E[i].patch], M));
        }
        y[i] = h * c;
    }

/* 2. berechne Integrale der goeberen Level aus denen der feineren */
    for (i = ne - nw - 1; i >= 0; i--) {
        y[i] = 0.5 * (y[E[i].son[0]] + y[E[i].son[1]] + y[E[i].son[2]] + y[E[i].son[3]]);
    }

/* 3. setze Integrale (f,psi) zusammen */
    for (i = 0; i < nw; i++) {
        c = 0;
        for (j = 0; j < W[i].element_number; j++)
            c += y[W[i].element[j]] * W[i].weight[j];
        (*rhs)[i] = c;
    }

/* Speicherplatz wieder freigeben */
    free_Gauss_Square(&Q, g + 1);
    free(y);
    return;
}

void WEMRHS2M(rhs, W, E, T, p, M, potential, g)
/* testet die Dirichlet-Daten des gegebenen Potentials */
double **rhs;                   /* zu berechnende rechte Seite                */
wavelet *W;                     /* Waveletliste                               */
element *E;                     /* hierarchische Elementliste                 */
vector3 ****T;                  /* Oberflaecheninterpolation                  */
unsigned int p;                 /* Anzahl der Patches                         */
unsigned int M;                 /* Zahl der Level                             */
unsigned int g;                 /* Quadratur-Grad auf dem feinsten Level      */
double *potential;
{
    unsigned int N = 1 << M;    /* N*N Elemente pro Patch auf dem Level M     */
    unsigned int nw;            /* Anzahl der Wavelets                        */
    unsigned int ne;            /* Anzahl der Elemente                        */
    signed int i, j;            /* Laufindizes durch die Wavelet/Elementliste */
    double h;                   /* Schrittweite                               */
    cubature *Q;                /* Kubatur-Formeln                            */
    double c;                   /* Wert der Integrals                         */
    vector2 t;                  /* Auswertepunkt der Gauss-Quadratur          */
    unsigned int k;             /* Laufindex fuer die Quadratur               */
    double *y;                  /* Array mit den Integralen (f,phi_m)         */
    unsigned int index;

    nw = p * N * N;             /* Anzahl der Wavelets */
    ne = p * (4 * N * N - 1) / 3;       /* Anzahl der Elemente */
    init_Gauss_Square(&Q, g + 1);       /* Kubatur-Formeln     */
    y = (double *) malloc(ne * sizeof(double));
    (*rhs) = (double *) malloc(nw * sizeof(double));

/* 1. Quadratur auf dem feinsten Level */
    h = 1. / N;
    for (i = ne - nw; i < ne; i++) {
        c = 0;
        for (k = 0; k < Q[g].nop; k++) {
            t.x = h * (E[i].index_s + Q[g].xi[k].x);
            t.y = h * (E[i].index_t + Q[g].xi[k].y);
            /*
             * Let's hope that these are the same as calculated beforehand!!!
             * c += Q[g].w[k] * f(Chi(t,T[E[i].patch],M));
             */
            index = (E[i].patch * N * N + E[i].index_t * N + E[i].index_s) * Q[g].nop + k;
            c += Q[g].w[k] * potential[index];
            //printf("%d\n", index);
        }
        y[i] = h * c;
    }

/* 2. berechne Integrale der goeberen Level aus denen der feineren */
    for (i = ne - nw - 1; i >= 0; i--) {
        y[i] = 0.5 * (y[E[i].son[0]] + y[E[i].son[1]] + y[E[i].son[2]] + y[E[i].son[3]]);
    }

/* 3. setze Integrale (f,psi) zusammen */
    for (i = 0; i < nw; i++) {
        c = 0;
        for (j = 0; j < W[i].element_number; j++)
            c += y[W[i].element[j]] * W[i].weight[j];
        (*rhs)[i] = c;
    }

/* Speicherplatz wieder freigeben */
    free_Gauss_Square(&Q, g + 1);
    free(y);
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

