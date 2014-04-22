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
 *  dwt.c  *
 ***********/


/*===========================================*
 *  Dieses Modul enthaelt alle Routinen der  *
 *  schnellen Wavelet-Transformationen.      *
 *===========================================*/


#include <stdlib.h>
#include <string.h>
#include "sparse.h"
#include "mask_pwl.h"
#include "dwt_pwl.h"


extern const unsigned int minLevel_pwl;


void multiple(C, Z, F, M, p, np)
/* Hilfsfunktion zur Wavelet-Transformation: Berechnet aus der 
   Elementliste F eine Liste der lokalen Gitterpunkte (Basisliste) 
   und eine Vielfachheitenliste der Punkte */
unsigned int ****C;             /* lokale Basisliste            */
unsigned int **Z;               /* Vielfachheitenliste          */
unsigned int **F;               /* Elementliste                 */
unsigned int M;                 /* (2^M*2^M) Elemente pro Patch */
unsigned int p;                 /* Anzahl der Patches           */
unsigned int np;                /* number of basic functions    */
{
    unsigned int N = 1 << M;    /* (N*N) Elemente pro Patch     */
    unsigned int zf;            /* Laufindex fuer Elelemtliste  */
    unsigned int zc;            /* Laufindex fuer Basisliste    */
    unsigned int i1, i2, i3;    /* Laufindizes fuer Punkte      */
    unsigned int arg;           /* temporaere Groesse           */

/* Speicherplatz allokieren */
    (*C) = (unsigned int ***) malloc(p * sizeof(unsigned int **));
    (*Z) = (unsigned int *) calloc(np, sizeof(unsigned int));

/* Vielfachheitenliste initialisieren */
    for (i1 = 0; i1 < p * N * N; i1++) {
        for (i2 = 0; i2 < 4; i2++)
            (*Z)[F[i1][i2]]++;
    }

/* durch feinstes Level gehen */
    zc = zf = 0;
    for (i1 = 0; i1 < p; i1++) {
        (*C)[i1] = (unsigned int **) malloc((N + 1) * sizeof(unsigned int *));
        for (i2 = 0; i2 <= N; i2++) {
            (*C)[i1][i2] = (unsigned int *) malloc((N + 1) * sizeof(unsigned int));
            for (i3 = 0; i3 <= N; i3++) {
                if ((i2 <= N / 2) && (i3 <= N / 2))
                    arg = F[zf++][0];
                else if ((i2 <= N / 2) && (i3 > N / 2))
                    arg = F[zf++][1];
                else if ((i2 > N / 2) && (i3 > N / 2))
                    arg = F[zf++][2];
                else
                    arg = F[zf++][3];

                (*C)[i1][i2][i3] = arg;
                if ((i2 != 0) && (i2 != N) && (i3 != 0) && (i3 != N)) {
                    (*Z)[arg] = 1;      /* ==> innerer Punkt */
                } else if (((i2 == 0) || (i2 == N)) ^ ((i3 == 0) || (i3 == N))) {
                    (*Z)[arg]--;        /* ==> Punkt auf Kante */
                } else;         /* ==> Punkt auf Ecke --> do nothing */

                if (i3 == N / 2)
                    zf--;
            }
            if (i2 == N / 2)
                zf -= N;
        }
    }
    return;
}


void dwtLin(a, F, M, p, np)
double *a;                      /* gesuchter/gegebener Koeffizientenvektor */
unsigned int **F;               /* Elementliste                            */
unsigned int M;                 /* (2^M*2^M) Elemente pro Patch            */
unsigned int p;                 /* Anzahl der Patches                      */
unsigned int np;                /* number of basic functions               */
{
    unsigned int ***C;          /* lokale Basisliste                       */
    unsigned int *Z;            /* Vielfachheitenliste                     */
    unsigned int m;             /* Laufindex fuer Level                    */
    unsigned int N;             /* (N*N) Elemente pro Patch auf Level M    */
    unsigned int n;             /* (n*n) Elemente pro Patch auf Level n    */
    unsigned int S;             /* Schrittweite zum naechsten Patch        */
    unsigned int i1, i2, i3;    /* Laufindizes durch die Basisliste        */
    sparse T, L;                /* Maskenmatrizen                          */
    unsigned int s, t;          /* Laufindizes durch die Maskenmatrizen    */
    double *b;                  /* Hilfsvektor                             */
    unsigned int arg;           /* temporaere Groesse                      */

/* SCHLEIFE UEBER DIE GITTER */
    N = 1 << M;
    multiple(&C, &Z, F, M, p, np);
    b = (double *) calloc(np, sizeof(double));
    for (m = M; m >= minLevel_pwl; m--) {
        dwt_mask_pwl(&T, &L, m, M);     /* berechne die Masken T und L      */
        n = 1 << m;             /* p*n*n Elemente auf Level m       */
        S = 1 << (M - m);       /* Schrittweite zum naechsten Punkt */

        /* 1. bilde Skalierungsfunktionen und Grund-Wavelets */
        for (i1 = 0; i1 < p; i1++) {
            for (i2 = 0; i2 <= n; i2++) {
                for (i3 = 0; i3 <= n; i3++) {
                    if (i3 % 2 == 1) {  /* Tensorprodukte psi_{m}(s)*phi_{m+1}(t) */
                        for (s = 0; s < T.row_number[i3]; s++) {
                            arg = C[i1][S * i2][S * T.index[i3][s]];
                            b[C[i1][S * i2][S * i3]] += T.value[i3][s] * a[arg];
                        }
                    } else {    /* Tensorprodukte phi_{m}(s)*phi_{m}(t) und phi_{m}(s)*psi_{m}(t) */
                        for (t = 0; t < T.row_number[i2]; t++) {
                            for (s = 0; s < T.row_number[i3]; s++) {
                                arg = C[i1][S * T.index[i2][t]][S * T.index[i3][s]];
                                b[C[i1][S * i2][S * i3]] += T.value[i3][s] * T.value[i2][t] * a[arg];
                            }
                        }
                    }
                }
            }
        }

        /* 2. a updaten */
        for (i1 = 0; i1 < p; i1++) {
            for (i2 = 0; i2 <= n; i2++) {
                for (i3 = 0; i3 <= n; i3++) {
                    a[C[i1][S * i2][S * i3]] = 0.5 * b[C[i1][S * i2][S * i3]];
                }
            }
        }

        /* 3. setze b = 0 im Falle einer Skalierungsfunktion */
        for (i1 = 0; i1 < p; i1++) {
            for (i2 = 0; i2 <= n; i2 += 2) {
                for (i3 = 0; i3 <= n; i3 += 2) {
                    b[C[i1][S * i2][S * i3]] = 0;
                }
            }
        }
        free_sparse(&T);
    }

/* Speicherplatz wieder freigeben */
    for (i1 = 0; i1 < p; i1++) {
        for (i2 = 0; i2 <= N; i2++)
            free(C[i1][i2]);
        free(C[i1]);
    }
    free(C);
    free(Z);
    free(b);
    return;
}


void tdwtLin(a, F, M, p, np)
double *a;                      /* gesuchter/gegebener Koeffizientenvektor */
unsigned int **F;               /* Elementliste                            */
unsigned int M;                 /* (2^M*2^M) Elemente pro Patch            */
unsigned int p;                 /* Anzahl der Patches                      */
unsigned int np;                /* number of basic functions               */
{
    unsigned int ***C;          /* lokale Basisliste                       */
    unsigned int *Z;            /* Punktvielfachheitenliste                */
    unsigned int m;             /* Laufindex fuer Level                    */
    unsigned int N;             /* (N*N) Elemente pro Patch auf Level M    */
    unsigned int n;             /* (n*n) Elemente pro Patch auf Level n    */
    signed int S;               /* Schrittweite zum naechsten Patch        */
    unsigned int i1, i2, i3;    /* Laufindizes durch die Basisliste        */
    sparse T, L;                /* Maskenmatrizen                          */
    unsigned int s, t;          /* Laufindizes durch die Maskenmatrizen    */
    double *b;                  /* Hilfsvektor                             */
    unsigned int arg;           /* temporaere Groesse                      */

/* SCHLEIFE UEBER DIE GITTER */
    N = 1 << M;
    multiple(&C, &Z, F, M, p, np);
    b = (double *) malloc(np * sizeof(double));
    for (m = minLevel_pwl; m <= M; m++) {
        dwt_mask_pwl(&T, &L, m, M);     /* berechne die Masken T und L      */
        n = 1 << m;             /* p*n*n Elemente auf Level m       */
        S = 1 << (M - m);       /* Schrittweite zum naechsten Punkt */

        /* 1. setze b = 0 im Falle einer Skalierungsfunktion */
        for (i1 = 0; i1 < p; i1++) {
            for (i2 = 0; i2 <= n; i2++) {
                for (i3 = 0; i3 <= n; i3++) {
                    b[C[i1][S * i2][S * i3]] = 0;
                }
            }
        }

        /* 2. bilde Skalierungsfunktionen */
        for (i1 = 0; i1 < p; i1++) {
            for (i2 = 0; i2 <= n; i2++) {
                for (i3 = 0; i3 <= n; i3++) {
                    if (i3 % 2 == 1) {  /* Tensorprodukte psi_{m}(s)*phi_{m+1}(t) */
                        for (s = 0; s < T.row_number[i3]; s++) {
                            arg = C[i1][S * i2][S * T.index[i3][s]];
                            b[arg] += T.value[i3][s] * a[C[i1][S * i2][S * i3]];
                        }
                    } else {    /* Tensorprodukte phi_{m}(s)*phi_{m}(t) und phi_{m}(s)*psi_{m}(t) */
                        for (t = 0; t < T.row_number[i2]; t++) {
                            for (s = 0; s < T.row_number[i3]; s++) {
                                arg = C[i1][S * T.index[i2][t]][S * T.index[i3][s]];
                                b[arg] += T.value[i3][s] * T.value[i2][t] * a[C[i1][S * i2][S * i3]];
                            }
                        }
                    }
                }
            }
        }

        /* 3. a updaten */
        for (i1 = 0; i1 < p; i1++) {
            for (i2 = 0; i2 <= n; i2++) {
                for (i3 = 0; i3 <= n; i3++) {
                    a[C[i1][S * i2][S * i3]] = 0.5 * b[C[i1][S * i2][S * i3]];
                }
            }
        }
        free_sparse(&T);
    }

/* Speicherplatz wieder freigeben */
    for (i1 = 0; i1 < p; i1++) {
        for (i2 = 0; i2 <= N; i2++)
            free(C[i1][i2]);
        free(C[i1]);
    }
    free(C);
    free(Z);
    free(b);
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

