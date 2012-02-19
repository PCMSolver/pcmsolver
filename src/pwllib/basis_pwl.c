/**************
 *  Basis2.c  *
 **************/


/*===================================================*
 *  Erzeugt die fuer das Wavelet-Galerkin-Verfahren  *
 *  notwendige Liste der Elemente ueber alle Gitter  *
 *  und die Waveletliste.                            *
 *===================================================*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "intvector_pwl.h"
#include "vector2.h"
#include "vector3.h"
#include "sparse.h"
#include "mask_pwl.h"
#include "kern.h"
#include "basis_pwl.h"


/*=================================================================*
 *  D E K L A R A T I O N   D E R   H I L F S F U N K T I O N E N  *
 *=================================================================*/

void generate_topology_pwl(unsigned int ****C, element_pwl * E, unsigned int p, unsigned int M, unsigned int nw);
/* Hilfsfunktion zur Wavelet-Transformation: Berechnet aus der 
   Elementliste E eine Liste der lokalen Gitterpunkte (Basisliste) */


void generate_canonical_single_scale_basis(wavelet_pwl * W, unsigned int ***C, element_pwl * E, unsigned int p, unsigned int m, unsigned int M);
/* erstellt die kanonische Einskalenbasis des Levels m */


void add_element(wavelet_pwl * w, element_pwl * E, double weight0, double weight1, double weight2, double weight3, unsigned int index);
/* fuegt zum Wavelet w das mit weight gewichtete Element index hinzu */


void add_wavelet(wavelet_pwl * w1, wavelet_pwl * w2, element_pwl * E, double weight);
/* fuegt zum Wavelet w1 das mit weight gewichtete Wavelet w2 hinzu */


/*===============================================================*
 *  D E F I N I T I O N   D E R   H I L F S F U N K T I O N E N  *
 *===============================================================*/

void unify_pwl(d, r, d1, r1, d2, r2)
/* bildet die Vereinigung K(d,r) = K(d1,r1) \cup K(d2,r2) */
vector3 *d, d1, d2;
double *r, r1, r2;
{
    vector3 z;
    double norm;

    z.x = d1.x - d2.x;
    z.y = d1.y - d2.y;
    z.z = d1.z - d2.z;
    norm = sqrt(z.x * z.x + z.y * z.y + z.z * z.z);

    if (norm + r2 <= r1) {      /* K(d2,r2) \subset K(d1,r1) */
        *d = d1;
        *r = r1;
    } else if (norm + r1 <= r2) {       /* K(d1,r1) \subset K(d2,r2) */
        *d = d2;
        *r = r2;
    } else {                    /* die Vereinigungsmenge ist keine Kugel */
        (*d).x = 0.5 * (d1.x + d2.x + (r1 - r2) / norm * z.x);
        (*d).y = 0.5 * (d1.y + d2.y + (r1 - r2) / norm * z.y);
        (*d).z = 0.5 * (d1.z + d2.z + (r1 - r2) / norm * z.z);
        *r = 0.5 * (r1 + r2 + norm);
    }
    return;
}


void generate_topology_pwl(C, E, p, M, nw)
/* Hilfsfunktion zur Wavelet-Transformation: Berechnet aus der 
   Elementliste E eine Liste der lokalen Gitterpunkte (Basisliste) */
unsigned int ****C;             /* lokale Basisliste            */
element_pwl *E;                 /* hierarchische Elementliste   */
unsigned int p;                 /* Anzahl der Patches           */
unsigned int M;                 /* (2^M*2^M) Elemente pro Patch */
unsigned int nw;                /* number of basic functions    */
{
    unsigned int N = 1 << M;    /* (N*N) Elemente pro Patch     */
    unsigned int ze;            /* Laufindex durch Elementliste */
    unsigned int i1, i2, i3;    /* Laufindizes fuer Punkte      */

/* Speicherplatz allokieren */
    (*C) = (unsigned int ***) malloc(p * sizeof(unsigned int **));

/* durch feinstes Level gehen */
    ze = p * (N * N - 1) / 3;
    for (i1 = 0; i1 < p; i1++) {
        (*C)[i1] = (unsigned int **) malloc((N + 1) * sizeof(unsigned int *));
        for (i2 = 0; i2 <= N; i2++) {
            (*C)[i1][i2] = (unsigned int *) malloc((N + 1) * sizeof(unsigned int));
            for (i3 = 0; i3 <= N; i3++) {
                if ((i2 <= N / 2) && (i3 <= N / 2))
                    (*C)[i1][i2][i3] = E[ze].vertex[0];
                else if ((i2 <= N / 2) && (i3 > N / 2))
                    (*C)[i1][i2][i3] = E[ze].vertex[1];
                else if ((i2 > N / 2) && (i3 > N / 2))
                    (*C)[i1][i2][i3] = E[ze].vertex[2];
                else
                    (*C)[i1][i2][i3] = E[ze].vertex[3];

                if (i3 != N / 2)
                    ze++;
            }
            if (i2 == N / 2)
                ze -= N;
        }
    }
    return;
}


void generate_canonical_single_scale_basis(W, C, E, p, m, M)
/* erstellt die kanonische Einskalenbasis des Levels m */
wavelet_pwl *W;                 /* Liste der kanonischen Einskalenbasis   */
unsigned int ***C;              /* lokale Basisliste                      */
element_pwl *E;                 /* hierarchische Elementliste             */
unsigned int p;                 /* Anzahl der Patches                     */
unsigned int m;                 /* 2^m*2^m Elemente pro Patch auf Level m */
unsigned int M;                 /* 2^M*2^M Elemente pro Patch auf Level M */
{
    unsigned int n = 1 << m;    /* p*n*n Elemente auf Level m             */
    unsigned int S = 1 << (M - m);      /* Schrittweite zum naechsten Patch       */
    unsigned int i1, i2, i3;    /* Laufindizes durch Gitter m             */
    unsigned int ze;            /* Zaehler fuer Elementliste              */

/* Initialisierung */
    for (i1 = 0; i1 < p; i1++) {
        for (i2 = 0; i2 <= n; i2++) {
            for (i3 = 0; i3 <= n; i3++) {
                W[C[i1][S * i2][S * i3]].element_number = 0;
                W[C[i1][S * i2][S * i3]].level = m;
            }
        }
    }

/* berechne die kanonische Einskalenbasis des Levels m */
    ze = p * (n * n - 1) / 3;   /* Index von Element 0 des Levels m */
    for (i1 = 0; i1 < p; i1++) {
        for (i2 = 0; i2 < n; i2++) {
            for (i3 = 0; i3 < n; i3++) {
                add_element(&W[C[i1][S * i2][S * i3]], E, 1, 0, 0, 0, ze);
                add_element(&W[C[i1][S * i2][S * (i3 + 1)]], E, 0, 1, 0, 0, ze);
                add_element(&W[C[i1][S * (i2 + 1)][S * (i3 + 1)]], E, 0, 0, 1, 0, ze);
                add_element(&W[C[i1][S * (i2 + 1)][S * i3]], E, 0, 0, 0, 1, ze);
                ze++;
            }
        }
    }
    return;
}


void add_element(w, E, weight0, weight1, weight2, weight3, index)
/* fuegt zum Wavelet w das mit weight gewichtete Element index hinzu */
wavelet_pwl *w;                 /* gegebenes Wavelet                        */
element_pwl *E;                 /* hierarchische Elementliste               */
double weight0;                 /* Gewicht im Eckpunkt 0 des Elements index */
double weight1;                 /* Gewicht im Eckpunkt 1 des Elements index */
double weight2;                 /* Gewicht im Eckpunkt 2 des Elements index */
double weight3;                 /* Gewicht im Eckpunkt 3 des Elements index */
unsigned int index;             /* Index des zu addierenden Elements        */
{
    unsigned int k;             /* Laufindex zur Suche nach Element index   */

    if ((weight0 == 0) && (weight1 == 0) && (weight2 == 0) && (weight3 == 0))
        return;                 /* es ist nix zu tun */

/* suche nach Element */
    for (k = 0; (k < w->element_number) && (w->element[k] != index); k++);
    if (k < w->element_number) {        /* Element schon vorhanden */
        w->weight[k][0] += weight0;
        w->weight[k][1] += weight1;
        w->weight[k][2] += weight2;
        w->weight[k][3] += weight3;
    } else {                    /* Element noch nicht vorhanden */
        if (w->element_number % delta == 0) {
            w->element = (unsigned int *) realloc(w->element, (w->element_number + delta) * sizeof(unsigned int));
            w->weight = (double (*)[4]) realloc(w->weight, (w->element_number + delta) * sizeof(double[4]));
        }

        w->element[w->element_number] = index;
        w->weight[w->element_number][0] = weight0;
        w->weight[w->element_number][1] = weight1;
        w->weight[w->element_number][2] = weight2;
        w->weight[w->element_number][3] = weight3;
        w->element_number++;
    }
    return;
}


void add_wavelet(w1, w2, E, weight)
/* fuegt zum Wavelet p1 das mit weight gewichtete Wavelet p2 hinzu */
wavelet_pwl *w1, *w2;           /* gegebene Wavelets                                */
element_pwl *E;                 /* hierarchische Elementliste                       */
double weight;                  /* Gewicht der zu addierenden Skalierungsfunktion   */
{
    unsigned int k;             /* Laufindizes zur Suche nach gemeinsamen Elementen */

    for (k = 0; k < w2->element_number; k++) {
        add_element(w1, E, weight * w2->weight[k][0], weight * w2->weight[k][1], weight * w2->weight[k][2], weight * w2->weight[k][3], w2->element[k]);
    }
    return;
}


/*===========================*
 *  E L E M E N T L I S T E  *
 *===========================*/

unsigned int generate_elementlist_pwl(E, P, F, p, M)
/* erstellt die hierarchsische Elementliste E */
element_pwl **E;                /* hierarchische Elementliste             */
vector3 *P;                     /* Punkteliste der Einskalenbasis         */
unsigned int **F;               /* Elementliste der Einskalenbasis        */
unsigned int p;                 /* Anzahl der Patches                     */
unsigned int M;                 /* 2^M*2^M Elemente pro Patch             */
{
    signed int m;               /* Laufindex fuer das Level               */
    unsigned int n;             /* p*n*n Elemente auf Level m             */
    unsigned int N = 1 << M;    /* p*N*N Elemente auf Level M             */
    unsigned int i1, i2, i3;    /* Laufindizes durch Gitter m             */
    unsigned int zf, ze;        /* Zaehler fuer Elementlisten             */
    unsigned int ne;            /* Laenge von E                           */
    vector3 d1, d2;             /* Umkreismittelpunkt von Element 0 und 2 */
    double r1, r2;              /* Umkreisradius      von Element 1 und 3 */

    ne = p * (4 * N * N - 1) / 3;       /* Laenge von E */
    (*E) = (element_pwl *) calloc(ne, sizeof(element_pwl));

/* initialisiere feinstes Gitter */
    zf = 0;
    ze = p * (N * N - 1) / 3;
    for (i1 = 0; i1 < p; i1++) {
        for (i2 = 0; i2 < N; i2++) {
            for (i3 = 0; i3 < N; i3++) {
                (*E)[ze].level = M;     /* Level des Elements        */
                (*E)[ze].patch = i1;    /* Patch des Elements        */
                (*E)[ze].index_s = i3;  /* Index in s-Richtung       */
                (*E)[ze].index_t = i2;  /* Index in t-Richtung       */
                (*E)[ze].vertex[0] = F[zf][0];  /* 1. Eckpunkt des Elements  */
                (*E)[ze].vertex[1] = F[zf][1];  /* 2. Eckpunkt des Elements  */
                (*E)[ze].vertex[2] = F[zf][2];  /* 3. Eckpunkt des Elements  */
                (*E)[ze].vertex[3] = F[zf][3];  /* 4. Eckpunkt des Elements  */

                /* setze die Soehne des Vaterelements */
                (*E)[ze].son[0] = (*E)[ze].son[1] = (*E)[ze].son[2] = (*E)[ze].son[3] = ne;

                /* bestimme Mittelpunkt und Radius des Elementumkreises */
                unify_pwl(&d1, &r1, P[F[zf][0]], 0, P[F[zf][2]], 0);
                unify_pwl(&d2, &r2, P[F[zf][1]], 0, P[F[zf][3]], 0);
                unify_pwl(&(*E)[ze].midpoint, &(*E)[ze].radius, d1, r1, d2, r2);

                ze++;
                zf++;
            }
        }
    }

/* Schleife ueber die groeberen Gitter */
    ze = ne;
    for (m = M - 1; m >= 0; m--) {
        n = 1 << m;
        ze -= 5 * p * n * n;
        zf = ze + p * n * n;
        for (i1 = 0; i1 < p; i1++) {
            for (i2 = 0; i2 < n; i2++) {
                for (i3 = 0; i3 < n; i3++) {
                    (*E)[ze].level = m; /* Level des Elements       */
                    (*E)[ze].patch = i1;        /* Patch des Elements       */
                    (*E)[ze].index_s = i3;      /* Index in s-Richtung      */
                    (*E)[ze].index_t = i2;      /* Index in t-Richtung      */
                    (*E)[ze].vertex[0] = (*E)[zf].vertex[0];    /* 1. Eckpunkt des Elements */
                    (*E)[ze].vertex[1] = (*E)[zf + 1].vertex[1];        /* 2. Eckpunkt des Elements */
                    (*E)[ze].vertex[2] = (*E)[zf + 2 * n + 1].vertex[2];        /* 3. Eckpunkt des Elements */
                    (*E)[ze].vertex[3] = (*E)[zf + 2 * n].vertex[3];    /* 4. Eckpunkt des Elements */
                    (*E)[ze].son[0] = zf;       /* 1. Sohn des Elements     */
                    (*E)[ze].son[1] = zf + 1;   /* 2. Sohn des Elements     */
                    (*E)[ze].son[2] = zf + 2 * n + 1;   /* 3. Sohn des Elements     */
                    (*E)[ze].son[3] = zf + 2 * n;       /* 4. Sohn des Elements     */

                    /* setze in den Soehnen das Vaterelement */
                    (*E)[zf].father = (*E)[zf + 1].father = (*E)[zf + 2 * n].father = (*E)[zf + 2 * n + 1].father = ze;

                    /* bestimme Mittelpunkt und Radius des Elementumkreises */
                    unify_pwl(&d1, &r1, (*E)[zf].midpoint, (*E)[zf].radius, (*E)[zf + 2 * n + 1].midpoint, (*E)[zf + 2 * n + 1].radius);
                    unify_pwl(&d2, &r2, (*E)[zf + 1].midpoint, (*E)[zf + 1].radius, (*E)[zf + 2 * n].midpoint, (*E)[zf + 2 * n].radius);
                    unify_pwl(&(*E)[ze].midpoint, &(*E)[ze].radius, d1, r1, d2, r2);

                    zf += 2;
                    ze++;
                }
                zf += 2 * n;
            }
        }
    }
    return (ne);
}


void complete_elementlist_pwl(W, E, p, M, nw)
/* erstellt die hierarchsische Elementliste E */
wavelet_pwl *W;                 /* Liste der Wavelets                 */
element_pwl *E;                 /* hierarchische Elementliste         */
unsigned int p;                 /* Anzahl der Patches                 */
unsigned int M;                 /* 2^M*2^M Elemente pro Patch         */
unsigned int nw;                /* Anzahl der Wavelets                */
{
    unsigned int ne;            /* Laenge hierarchische Elementliste  */
    unsigned int i, j;          /* Laufindizes durch die Waveletliste */
    unsigned int k;             /* zu bearbeitendes Element           */

/* Elementliste updaten */
    for (i = 0; i < nw; i++) {
        for (j = 0; j < W[i].element_number; j++) {
            k = W[i].element[j];
            if (E[k].wavelet_number % delta == 0) {
                E[k].wavelet = (unsigned int *) realloc(E[k].wavelet, (E[k].wavelet_number + delta) * sizeof(unsigned int));
            }
            E[k].wavelet[E[k].wavelet_number++] = i;
        }
    }

/* Speicherplatz der Elementliste optimieren */
    ne = p * (4 * (1 << 2 * M) - 1) / 3;        /* Laenge von E */
    for (i = 0; i < ne; i++)
        E[i].wavelet = (unsigned int *) realloc(E[i].wavelet, E[i].wavelet_number * sizeof(unsigned int));
    return;
}


void free_elementlist_pwl(E, p, M)
/* gibt den Speicherplatz der hierarchischen Elementliste E frei */
element_pwl **E;                /* hierarchische Elementliste */
unsigned int p;                 /* Anzahl der Patches         */
unsigned int M;                 /* 2^M*2^M Elemente pro Patch */
{
    unsigned int ne;            /* Laenge von E               */
    unsigned int i;             /* Laufindex durch E          */

    ne = p * (4 * (1 << 2 * M) - 1) / 3;        /* Anzahl aller Elemente */

    for (i = 0; i < ne; i++)
        free((*E)[i].wavelet);
    free(*E);
}


/*===========================*
 *  W A V E L E T L I S T E  *
 *===========================*/

void generate_waveletlist_pwl(W, E, p, M, nw)
/* erstellt die Waveletliste W */
wavelet_pwl **W;                /* Liste der Wavelets                   */
element_pwl *E;                 /* hierarchische Elementliste           */
unsigned int p;                 /* Anzahl der Patches                   */
unsigned int M;                 /* 2^M*2^M Elemente pro Patch           */
unsigned int nw;                /* Laenge von W                         */
{
    unsigned int ***C;          /* lokale Basisliste                    */
    unsigned int N = 1 << M;    /* p*N*N Elemente auf Level M           */
    unsigned int n;             /* p*n*n Elemente auf Level m           */
    signed int m;               /* Laufindex fuer das Level             */
    unsigned int i1, i2, i3;    /* Laufindizes durch Gitter m           */
    wavelet_pwl *w;             /* Zeiger auf zu bearbeitendes Wavelet  */
    sparse T, L;                /* Maskenmatrizen                       */
    unsigned int s, t;          /* Laufindizes durch die Maskenmatrizen */
    unsigned int S;             /* Schrittweite zum naechsten Patch     */
    unsigned int arg;           /* temporaere Groesse                   */
    wavelet_pwl *G;             /* temporaere Waveletliste              */

/* 1. Initialisierung */
    generate_topology_pwl(&C, E, p, M, nw);
    G = (wavelet_pwl *) calloc(nw, sizeof(wavelet_pwl));
    (*W) = (wavelet_pwl *) calloc(nw, sizeof(wavelet_pwl));

/* 2. Schleife ueber die Gitter */
    for (m = M; m >= minLevel_pwl; m--) {
        dwt_mask_pwl(&T, &L, m, m);     /* berechne die Masken T und L      */
        n = 1 << m;             /* p*n*n Elemente auf Level m       */
        S = 1 << (M - m);       /* Schrittweite zum naechsten Punkt */

        /* 2.1. bilde Wavelets */
        generate_canonical_single_scale_basis(G, C, E, p, m, M);
        for (i1 = 0; i1 < p; i1++) {
            for (i2 = 0; i2 <= n; i2++) {
                for (i3 = 0; i3 <= n; i3++) {
                    w = &(*W)[C[i1][S * i2][S * i3]];   /* zu bearbeitendes Wavelet */
                    w->level = m;       /* Level des Wavelets       */

                    /* fuege Elemente des Wavelets hinzu */
                    if (i3 % 2 == 1) {  /* Tensorprodukte psi_{m}(s)*phi_{m+1}(t) */
                        for (s = 0; s < T.row_number[i3]; s++) {
                            arg = C[i1][S * i2][S * T.index[i3][s]];
                            add_wavelet(w, &G[arg], E, 0.5 * T.value[i3][s]);
                        }
                    } else if (i2 % 2 == 1) {   /* Tensorprodukte phi_{m}(s)*psi_{m}(t) */
                        for (t = 0; t < T.row_number[i2]; t++) {
                            for (s = 0; s < T.row_number[i3]; s++) {
                                arg = C[i1][S * T.index[i2][t]][S * T.index[i3][s]];
                                add_wavelet(w, &G[arg], E, 0.5 * T.value[i3][s] * T.value[i2][t]);
                            }
                        }
                    }

                    /* bestimme Soehne des Wavelets */
                    if ((i3 % 4 == 2) && (i2 % 2 == 0)) {
                        if (i2 < n) {
                            w->son_number = 4;
                            w->son = (unsigned int *) realloc(w->son, 4 * sizeof(unsigned int));
                            w->son[0] = C[i1][S * i2][S * (i3 - 1)];
                            w->son[1] = C[i1][S * i2][S * (i3 + 1)];
                            w->son[2] = C[i1][S * (i2 + 1)][S * (i3 - 1)];
                            w->son[3] = C[i1][S * (i2 + 1)][S * (i3 + 1)];
                        } else {
                            w->son_number = 2;
                            w->son = (unsigned int *) realloc(w->son, 2 * sizeof(unsigned int));
                            w->son[0] = C[i1][S * i2][S * (i3 - 1)];
                            w->son[1] = C[i1][S * i2][S * (i3 + 1)];
                        }
                    } else if ((i2 % 4 == 2) && (i3 % 4 == 0)) {
                        if (i3 < n) {
                            w->son_number = 4;
                            w->son = (unsigned int *) realloc(w->son, 4 * sizeof(unsigned int));
                            w->son[0] = C[i1][S * (i2 - 1)][S * i3];
                            w->son[1] = C[i1][S * (i2 - 1)][S * (i3 + 2)];
                            w->son[2] = C[i1][S * (i2 + 1)][S * i3];
                            w->son[3] = C[i1][S * (i2 + 1)][S * (i3 + 2)];
                        } else {
                            w->son_number = 2;
                            w->son = (unsigned int *) realloc(w->son, 2 * sizeof(unsigned int));
                            w->son[0] = C[i1][S * (i2 - 1)][S * i3];
                            w->son[1] = C[i1][S * (i2 + 1)][S * i3];
                        }
                    }
                }
            }
        }
        free_sparse(&T);
    }

/* 3. Skalierungsfunktionen hinzufuegen */
    generate_canonical_single_scale_basis(*W, C, E, p, minLevel_pwl - 1, M);

/* 4. Speicherplatz wieder freigeben */
    for (i1 = 0; i1 < nw; i1++) {
        free(G[i1].element);
        free(G[i1].weight);
    }
    for (i1 = 0; i1 < p; i1++) {
        for (i2 = 0; i2 <= N; i2++)
            free(C[i1][i2]);
        free(C[i1]);
    }
    free(G);
    free(C);
    return;
}


void set_quadrature_level_pwl(W, E, p, M, nw)
/* verfeinert Grobgitterelemente */
wavelet_pwl *W;                 /* Liste der Wavelets                     */
element_pwl *E;                 /* hierarchische Elementliste             */
unsigned int p;                 /* Anzahl der Patches                     */
unsigned int M;                 /* 2^M*2^M Elemente pro Patch             */
unsigned int nw;                /* Laenge von W                           */
{
    unsigned int i, j, k, l;    /* Laufindizes durch Wavelet/Elementliste */
    unsigned int ind;           /* Index des zu untersuchenden Elements   */
    unsigned int minLevel_pwl;      /* Minimales Level fuer die Quadratur     */
    unsigned int noe;           /* number of elements eines Wavelets      */
    double w[4];                /* Gewichte der Elementfunktionen         */

    minLevel_pwl = min_quadrature_level;    /* minimales Quadraturlevel */
    if (minLevel_pwl > M)
        minLevel_pwl = M;           /* sonst gibt's Aerger      */

    for (i = 0; (i < nw) && (W[i].level <= minLevel_pwl); i++) {
        for (j = W[i].level; (j <= minLevel_pwl); j++) {
            noe = W[i].element_number;
            for (k = 0; k < noe; k++) {
                ind = W[i].element[k];  /* zu untersuchendes Element  */
                if (E[ind].level < minLevel_pwl) {  /* also zu grob -> verfeinere */
                    W[i].element[k] = E[ind].son[0];
                    memcpy(w, W[i].weight[k], 4 * sizeof(double));
                    W[i].element = (unsigned int *) realloc(W[i].element, (W[i].element_number + 3) * sizeof(unsigned int));
                    W[i].weight = (double (*)[4]) realloc(W[i].weight, (W[i].element_number + 3) * sizeof(double[4]));

                    /* neue Gewichte bestimmen */
                    for (l = 1; l < 4; l++) {
                        W[i].element[W[i].element_number] = E[ind].son[l];
                        W[i].weight[W[i].element_number][l - 1] = 0.25 * (w[l - 1] + w[l]);
                        W[i].weight[W[i].element_number][l] = 0.5 * w[l];
                        W[i].weight[W[i].element_number][(l + 1) % 4] = 0.25 * (w[l] + w[(l + 1) % 4]);
                        W[i].weight[W[i].element_number][(l + 2) % 4] = 0.125 * (w[0] + w[1] + w[2] + w[3]);
                        W[i].element_number++;
                    }

                    W[i].element[k] = E[ind].son[0];
                    W[i].weight[k][0] = 0.5 * w[0];
                    W[i].weight[k][1] = 0.25 * (w[0] + w[1]);
                    W[i].weight[k][2] = 0.125 * (w[0] + w[1] + w[2] + w[3]);
                    W[i].weight[k][3] = 0.25 * (w[3] + w[0]);
                }
            }
        }
    }
    return;
}


void simplify_waveletlist_pwl(W, E, p, M, nw)
/* optimiert die Waveletliste W */
wavelet_pwl *W;                 /* Liste der Wavelets                          */
element_pwl *E;                 /* hierarchische Elementliste                  */
unsigned int p;                 /* Anzahl der Patches                          */
unsigned int M;                 /* 2^M*2^M Elemente pro Patch                  */
unsigned int nw;                /* Laenge von W                                */
{
    unsigned int i, k, l;       /* Laufindizes durch die Wavelet/Elementliste  */
    signed int j;               /* Laufindex durch die Wavelet/Elementliste    */
    unsigned int s0, s1, s2, s3;        /* Laufindizes fuer die Suche nach den Soehnen */
    unsigned int noe;           /* number of elements eines Wavelets           */
    unsigned int *prototype;    /* Liste fuer den Zugriff Prototyp -> Wavelet  */
    unsigned int prototype_number;      /* Anzahl der diversen Protptypen              */
    unsigned int minLevel_pwl;      /* Minimales Level fuer die Quadratur          */

/* 1. Untersuche, ob Feingitterelemente durch ihren Vater ersetzt werden   */
/*    koennen. Falls ja, tue dies, und setze das Gewicht der Soehne auf 0. */
    minLevel_pwl = min_quadrature_level;    /* minimales Quadraturlevel */
    for (i = 0; i < nw; i++) {
        if (W[i].level > minLevel_pwl) {
            noe = W[i].element_number;  /* Zahl der Eintraege in der Elementliste des Wavelets */
            for (s0 = 0; s0 < noe; s0++) {
                j = W[i].element[s0];   /* zu untersuchendes Element */
                if (E[j].level == W[i].level) { /* also Feingitterelement -> suche nach allen Soehnen *//* ueberpruefe, ob zu untersuchendes Element der 0. Sohn ist */
                    j = E[j].father;    /* Vaterelement */
                    if (W[i].element[s0] == E[j].son[0]) {      /* 0. Sohn vorhanden -> ueberpruefe auf 1. Sohn */
                        for (s1 = 0; (s1 < noe) && (E[j].son[1] != W[i].element[s1]); s1++);
                        if ((s1 < noe) && (fabs(2 * W[i].weight[s0][1] - W[i].weight[s0][0] - W[i].weight[s1][1]) < eps)) {     /* 0.+1. Sohn vorhanden -> ueberpruefe auf 2. Sohn */
                            for (s2 = 0; (s2 < noe) && (E[j].son[2] != W[i].element[s2]); s2++);
                            if ((s2 < noe) && (fabs(2 * W[i].weight[s1][2] - W[i].weight[s1][1] - W[i].weight[s2][2]) < eps)) { /* 0.+1.+2. Sohn vorhanden -> ueberpruefe auf 3. Sohn */
                                for (s3 = 0; (s3 < noe) && (E[j].son[3] != W[i].element[s3]); s3++);
                                {
                                    if ((s3 < noe) && (fabs(2 * W[i].weight[s2][3] - W[i].weight[s2][2] - W[i].weight[s3][3]) < eps)
                                        && (fabs(2 * W[i].weight[s3][0] - W[i].weight[s3][3] - W[i].weight[s0][0]) < eps)
                                        && (fabs(4 * W[i].weight[s0][2] - W[i].weight[s0][0] - W[i].weight[s1][1] - W[i].weight[s2][2] - W[i].weight[s3][3]) < eps)) {  /* alle 4 Soehne sind vorhanden und haben gleiche Gewichte */
                                        W[i].element[s0] = j;
                                        W[i].weight[s0][0] = 2 * W[i].weight[s0][0];
                                        W[i].weight[s0][1] = 2 * W[i].weight[s1][1];
                                        W[i].weight[s0][2] = 2 * W[i].weight[s2][2];
                                        W[i].weight[s0][3] = 2 * W[i].weight[s3][3];
                                        memset(W[i].weight[s1], 0, 4 * sizeof(double));
                                        memset(W[i].weight[s2], 0, 4 * sizeof(double));
                                        memset(W[i].weight[s3], 0, 4 * sizeof(double));
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /* Groesse der Elementliste anpassen */
            k = 0;
            while (k < W[i].element_number) {
                if ((W[i].weight[k][0] == 0) && (W[i].weight[k][1] == 0) && (W[i].weight[k][2] == 0) && (W[i].weight[k][3] == 0)) {
                    W[i].element_number--;      /* streiche Element */
                    for (l = k; l < W[i].element_number; l++) {
                        W[i].element[l] = W[i].element[l + 1];
                        memcpy(W[i].weight[l], W[i].weight[l + 1], 4 * sizeof(double));
                    }
                } else
                    k++;
            }

            W[i].element = (unsigned int *) realloc(W[i].element, W[i].element_number * sizeof(unsigned int));
            W[i].weight = (double (*)[4]) realloc(W[i].weight, W[i].element_number * sizeof(double[4]));
        }
    }

/* 2. bestimme die Liste der Prototypen */
    prototype = NULL;
    prototype_number = 0;

    for (i = 0; i < nw; i++) {
        for (j = prototype_number - 1; j >= 0; j--) {
            if (W[i].element_number == W[prototype[j]].element_number) {
                /* gleiche Elementanzahl: gleiche Gewichte? */
                for (k = 0; (k < W[i].element_number) && (fabs(W[i].weight[k][0] - W[prototype[j]].weight[k][0]) < eps)
                     && (fabs(W[i].weight[k][1] - W[prototype[j]].weight[k][1]) < eps)
                     && (fabs(W[i].weight[k][2] - W[prototype[j]].weight[k][2]) < eps)
                     && (fabs(W[i].weight[k][3] - W[prototype[j]].weight[k][3]) < eps); k++);

                if (k == W[i].element_number) { /* alle Gewichte gleich ==> ersetze die Liste */
                    free(W[i].weight);  /* der Gewichte durch die des Prototypen und  */
                    W[i].weight = W[prototype[j]].weight;       /* und verlasse die Schleife bzgl. j.         */
                    break;
                }
            }
        }

        if (j == -1) {          /* keinen passenden Prototyp gefunden */
            if (prototype_number % delta == 0)
                prototype = (unsigned int *) realloc(prototype, (prototype_number + delta) * sizeof(unsigned int));
            prototype[prototype_number++] = i;
        }
    }

    free(prototype);
    printf("%d prototypes\n", prototype_number);
    return;
}


void free_waveletlist_pwl(W, nw)
/* gibt den Speicherplatz der Waveletliste W frei */
wavelet_pwl **W;                /* Liste der Wavelets  */
unsigned int nw;                /* Anzahl der Wavelets */
{
    unsigned int i;             /* Laufindex durch W                          */
    unsigned int *prototype;    /* Liste fuer den Zugriff Prototyp -> Wavelet */
    signed int prototype_number;        /* Anzahl der diversen Protptypen             */
    signed int j;               /* Laufindex durch die Liste der Prototypen   */

    prototype = NULL;
    prototype_number = 0;

/* bestimme die Liste der Prototypen */
    for (i = 0; i < nw; i++) {
        free((*W)[i].son);
        free((*W)[i].element);

        /* suche nach passendem Prototypen */
        for (j = prototype_number - 1; j >= 0; j--)
            if ((*W)[i].weight == (*W)[prototype[j]].weight)
                break;

        if (j == -1) {          /* keinen passenden Prototyp gefunden */
            if (prototype_number % delta == 0)
                prototype = (unsigned int *) realloc(prototype, (prototype_number + delta) * sizeof(unsigned int));
            prototype[prototype_number++] = i;
        }
    }

/* Speicherplatz der Prototypen freigeben */
    for (i = 0; i < prototype_number; i++)
        free((*W)[prototype[i]].weight);
    free(prototype);
    free(*W);
}


void print_waveletlist_pwl(W, nw)
/* gibt die in der Waveletliste W definierten Wavelets aus */
wavelet_pwl *W;                 /* Waveletliste */
unsigned int nw;                /* Laenge von W */
{
    unsigned int i, j;          /* Laufindizes  */

    for (i = 0; i < nw; i++) {
        printf("\nWavelet: %d \nLevel:   %d \n%d Soehne:   ", i, W[i].level, W[i].son_number);
        for (j = 0; j < W[i].son_number; j++)
            printf("%d ", W[i].son[j]);
        printf("\n%d Elemente: ", W[i].element_number);
        for (j = 0; j < W[i].element_number; j++) {
            printf("[%d|(%g,%g,%g,%g)] ", W[i].element[j], W[i].weight[j][0], W[i].weight[j][1], W[i].weight[j][2], W[i].weight[j][3]);
        }
        printf("\n");
    }

    return;
}


/*===================================*
 *  A B S T A N D S F U N K T I O N  *
 *===================================*/

double distance_pwl(element1, element2)
/* Berechnet den Abstand zwischen den Elementen element1 und element2 */
element_pwl *element1, *element2;       /* Pointer auf die zwei Elemente */
{
    double dx, dy, dz;          /* x/y/z-Abstand zweier Elemente */
    double c1, c2;              /* Skalierungsfaktoren           */
    double dist;                /* zu berechnender Abstand       */

    dx = element1->midpoint.x - element2->midpoint.x;
    dy = element1->midpoint.y - element2->midpoint.y;
    dz = element1->midpoint.z - element2->midpoint.z;
    dist = sqrt(dx * dx + dy * dy + dz * dz) - element1->radius - element2->radius;
    c1 = element1->radius * (1 << element1->level);
    c2 = element2->radius * (1 << element2->level);

    if (c1 < c2)
        return (scaling_factor * dist / c2);
    else
        return (scaling_factor * dist / c1);
}
