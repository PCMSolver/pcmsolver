/*******************
 *  Compression.c  *
 *******************/


/*===========================================*
 *  Bestimmt in der Steifigkeitsmatrix alle  *
 *  nach der 1. und 2. Kompression zu        *
 *  berechnenden Eintraege.                  *
 *===========================================*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "sparse2.h"
#include "intvector_pwl.h"
#include "vector2.h"
#include "vector3.h"
#include "basis_pwl.h"
#include "compression_pwl.h"


/*========================================================*
 *  Deklaration der Hilfsfunktionen und der bounding box  *
 *========================================================*/

typedef struct {                /* Typdeklaration bounding box */
    double r_x, r_y, r_z;
    double m_x, m_y, m_z;
} bounding_box;


void compute_bounding_boxes(bounding_box ** B, wavelet_pwl * W, element_pwl * E, unsigned int np);


unsigned int wavelet_wavelet_criterion_pwl(bounding_box * B, wavelet_pwl * W, element_pwl * E, unsigned int ind1, unsigned int ind2, double c1, double c2);


/*==================================================*
 *  Berechnet die bounding boxes fuer die Wavelets  *
 *==================================================*/

void compute_bounding_boxes(B, W, E, np)
bounding_box **B;               /* bounding boxes             */
wavelet_pwl *W;                 /* Liste der Wavelets         */
element_pwl *E;                 /* hierarchische Elementliste */
unsigned int np;                /* Anzahl der Basisfunktionen */
{
    unsigned int i, j, k;       /* Laufindizes                */
    double min_x, min_y, min_z; /* minimale Ausdehnung        */
    double max_x, max_y, max_z; /* maximale Ausdehnung        */

    (*B) = (bounding_box *) malloc(np * sizeof(bounding_box));
    for (i = 0; i < np; i++) {
        min_x = max_x = E[W[i].element[0]].midpoint.x;
        min_y = max_y = E[W[i].element[0]].midpoint.y;
        min_z = max_z = E[W[i].element[0]].midpoint.z;
        for (j = 0; j < W[i].element_number; j++) {
            k = W[i].element[j];
            if (E[k].midpoint.x - E[k].radius < min_x)
                min_x = E[k].midpoint.x - E[k].radius;
            if (E[k].midpoint.y - E[k].radius < min_y)
                min_y = E[k].midpoint.y - E[k].radius;
            if (E[k].midpoint.z - E[k].radius < min_z)
                min_z = E[k].midpoint.z - E[k].radius;
            if (E[k].midpoint.x + E[k].radius > max_x)
                max_x = E[k].midpoint.x + E[k].radius;
            if (E[k].midpoint.y + E[k].radius > max_y)
                max_y = E[k].midpoint.y + E[k].radius;
            if (E[k].midpoint.z + E[k].radius > max_z)
                max_z = E[k].midpoint.z + E[k].radius;
        }
        (*B)[i].r_x = 0.5 * (max_x - min_x);
        (*B)[i].r_y = 0.5 * (max_y - min_y);
        (*B)[i].r_z = 0.5 * (max_z - min_z);
        (*B)[i].m_x = 0.5 * (max_x + min_x);
        (*B)[i].m_y = 0.5 * (max_y + min_y);
        (*B)[i].m_z = 0.5 * (max_z + min_z);
    }
    return;
}


/*=========================================================*
 *  wavelet_wavelet_criterion_pwl: 1.+2. Kompression	   *
 *  result = 0: T(ind,ind2) = 0 gemaess 1./2. Kompression  *
 *  result = 1: sonst 					   *
 *=========================================================*/

unsigned int wavelet_wavelet_criterion_pwl(B, W, E, ind1, ind2, c1, c2)
bounding_box *B;                /* bounding boxes                                */
wavelet_pwl *W;                 /* Liste der Wavelets                            */
element_pwl *E;                 /* hierarchische Elementliste                    */
unsigned int ind1, ind2;        /* Indizes der beiden Wavelets                   */
double c1, c2;                  /* Abschneideparameter 1./2.  Kompression        */
{
    unsigned int i, j;          /* Laufindizes durch die Elemente                */
    unsigned int k, l;          /* Indizes der Elemente                          */
    double h1, h2;              /* halbe Kantenlaenge von Element1 bzw. Element2 */
    double s1, t1;              /* (s1,t1) = Mittelpunkt von Element1            */
    double s2, t2;              /* (s2,t2) = Mittelpunkt von Element2            */
    double dist;                /* Abstand zwischen den MP zweier Elemente       */
    double dx, dy, dz;          /* x/y/z-Abstand zweier Wavelets                 */

/* 1. Kompression: nehme die bounding boxes um die Wavelets */
    dx = fabs(B[ind1].m_x - B[ind2].m_x) - B[ind1].r_x - B[ind2].r_x;
    dy = fabs(B[ind1].m_y - B[ind2].m_y) - B[ind1].r_y - B[ind2].r_y;
    dz = fabs(B[ind1].m_z - B[ind2].m_z) - B[ind1].r_z - B[ind2].r_z;
    if (dx < 0)
        dx = 0;
    if (dy < 0)
        dy = 0;
    if (dz < 0)
        dz = 0;
    if (dx * dx + dy * dy + dz * dz >= c1 * c1)
        return (0);

/* 2. Kompression: teste einzelne Elemente */
    for (i = 0; i < W[ind1].element_number; i++) {
        k = W[ind1].element[i];
        h1 = 0.5 / (1 << E[k].level);
        s1 = h1 * (2 * E[k].index_s + 1);
        t1 = h1 * (2 * E[k].index_t + 1);

        for (j = 0; j < W[ind2].element_number; j++) {
            l = W[ind2].element[j];

            /* entweder Elemente liegen auf dem selben Patch */
            if (E[k].patch == E[l].patch) {
                h2 = 0.5 / (1 << E[l].level);
                s2 = h2 * (2 * E[l].index_s + 1);
                t2 = h2 * (2 * E[l].index_t + 1);

                /* berechne Unendlichnorm des Abstandes der beiden Mittelpunkte */
                dist = (fabs(s1 - s2) < fabs(t1 - t2)) ? fabs(t1 - t2) : fabs(s1 - s2);

                /* Abstand der Elemente */
                if (fabs(fabs(dist - h2) - h1) < c2)
                    return (1);
            }

            /* oder Elemente liegen auf verschiedenen Patches */
            else if (distance_pwl(&E[k], &E[l]) < c2)
                return (1);
        }
    }
    return (0);
}


/*=================*
 *  Hauptprogramm  *
 *=================*/

double compression_pwl(T, W, E, p, M, np)
sparse2 *T;                     /* komprimierte Steifigkeitsmatrix               */
wavelet_pwl *W;                 /* Liste der Wavelets                            */
element_pwl *E;                 /* hierarchische Elementliste                    */
unsigned int p;                 /* Zahl der Patches                              */
unsigned int M;                 /* 2^M*2^M Patches pro Parametergebiet           */
unsigned int np;                /* Anzahl der Basisfunktionen                    */
{
    bounding_box *B;            /* bounding boxes fuer die Wavelets              */
    unsigned int m1, m2;        /* Laufindizes durch Level                       */
    unsigned int ind1, ind2;    /* Argumente der Vater-Wavelets                  */
    unsigned int i, j;          /* Spaltenlaufindizes durch Steifigkeitsmatrix T */
    unsigned int k, l;          /* Laufindizes durch Soehne                      */
    unsigned int *rn;           /* Zaehlvektor zur Symmetrisierung               */
    double max_radius;          /* maximaler Umkreis der Elemente am Level 0     */
    double **c1, **c2;          /* Abschneideparameter in 1./2. Kompression      */
    double d1, d2;              /* Hilfskonstanten                               */
    unsigned int nnz;           /* nunber of nonzero elements of T               */

/* berechne Abschneideparameter */
    max_radius = 0;
    for (i = 0; E[i].level == 0; i++)
        if (max_radius < E[i].radius)
            max_radius = E[i].radius;
    c1 = (double **) malloc((M + 1) * sizeof(double *));
    c2 = (double **) malloc((M + 1) * sizeof(double *));
    for (i = 0; i <= M; i++) {
        c1[i] = (double *) malloc((i + 1) * sizeof(double));
        c2[i] = (double *) malloc((i + 1) * sizeof(double));
        for (j = 0; j <= i; j++) {
            c1[i][j] = a * pow(2, (M * (2 * dp - op) - (i + j) * (dp + td_pwl)) / (2 * td_pwl + op));
            d1 = pow(2, (M * (2 * dp - op) - (i + j) * dp - i * td_pwl) / (td_pwl + op));       /* alter Parameter */
            d2 = pow(2, (2 * M - (i + j)) * (2 * dp - op) / ((2 * td_pwl + op) * (td_pwl + op))) * pow(2, (M * (2 * dp - op) - i * (dp + td_pwl + 1) - j * (dp - 1)) / (td_pwl + op));
            if (d1 > d2)
                c2[i][j] = a * d1;
            else
                c2[i][j] = a * d2;
            if (c1[i][j] < a / (1 << j))
                c1[i][j] = a / (1 << j);
            if (c2[i][j] < a / (1 << i))
                c2[i][j] = a / (1 << i);
            if ((i < minLevel_pwl) || (j < minLevel))
                c1[i][j] = c2[i][j];
            c1[i][j] *= max_radius / scaling_factor;    /* Gebiet relativieren */
        }
    }

/* Initialisierung: belege den Block 
   (minLevel_pwl-1:minLevel,minLevel-1:minLevel) mit 1 */
debugFile = fopen("debug.out","a");
fprintf(debugFile,">>> COMP CC");
for(i  = 0; i < M+1; ++i){
  for(j = 0; j < i+1; ++j){
    fprintf(debugFile, "[%lf %lf]", c1[i][j], c2[i][j]);
  }
  fprintf(debugFile,"\n");
}
fprintf(debugFile,"\n<<< COMP CC\n");
fflush(debugFile);
fclose(debugFile);
    init_pattern(T, np, np, 20);
    compute_bounding_boxes(&B, W, E, np);
    for (i = 0; (i < np) && (W[i].level < minLevel_pwl); i++) {
        for (j = 0; j <= i; j++)
            set_pattern(T, i, j);
    }

    for (i = 0; W[i].level < M; i++) {
        m1 = W[i].level + 1;

        /* bestimme aus den Bloecken (m1-1,1:m1-1) die Bloecke (m1,1:m1) */
        for (j = 0; j < T->row_number[i]; j++) {
            ind2 = T->index[i][j];
            m2 = W[ind2].level;
            for (k = 0; k < W[i].son_number; k++) {
                ind1 = W[i].son[k];
                if (wavelet_wavelet_criterion_pwl(B, W, E, ind1, ind2, c1[m1][m2], c2[m1][m2])) {
                    set_pattern(T, ind1, ind2);
                    if (m1 == m2 + 1) {
                        for (l = 0; l < W[ind2].son_number; l++) {
                            if ((W[ind2].son[l] <= ind1) && (wavelet_wavelet_criterion_pwl(B, W, E, ind1, W[ind2].son[l], c1[m1][m1], c2[m1][m1]))) {
                                set_pattern(T, ind1, W[ind2].son[l]);
                            }
                        }
                    }
                }
            }

            /* wegen halbem Diagonalblock */
            if (m1 == m2 + 1) {
                for (k = 0; k < W[ind2].son_number; k++) {
                    if (wavelet_wavelet_criterion_pwl(B, W, E, W[ind2].son[k], i, c1[m1][m2], c2[m1][m2])) {
                        set_pattern(T, W[ind2].son[k], i);
                        {
                            for (l = 0; l < W[i].son_number; l++) {
                                if ((W[i].son[l] <= W[ind2].son[k]) && (wavelet_wavelet_criterion_pwl(B, W, E, W[ind2].son[k], W[i].son[l], c1[m1][m1], c2[m1][m1]))) {
                                    set_pattern(T, W[ind2].son[k], W[i].son[l]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

/* Symmetrisierung der Steifigkeitsmatrix */
    nnz = 0;
    rn = (unsigned int *) calloc(np, sizeof(unsigned int));
    for (i = 0; i < np; i++) {
        rn[i] += T->row_number[i];
        for (j = 0; j + 1 < T->row_number[i]; j++)
            rn[T->index[i][j]]++;
    }
    for (i = 0; i < np; i++) {
        T->index[i] = (unsigned int *) realloc(T->index[i], (rn[i] + 1) * sizeof(unsigned int));
        T->max_row_number[i] = rn[i];
        T->index[i][rn[i]] = np;        /* setze Waechterelement */
        nnz += rn[i];
    }
    for (i = 0; i < np; i++) {
        for (j = 0; j + 1 < T->row_number[i]; j++)
            T->index[T->index[i][j]][T->row_number[T->index[i][j]]++] = i;
    }

/* Speicherplatz wieder freigeben */
/*
printf("A-priori Compression:            %.5f %%\n",100.0*nnz/np/np);
*/
    finish_pattern(T);
    for (i = 0; i <= M; i++) {
        free(c1[i]);
        free(c2[i]);
    }
    free(c1);
    free(c2);
    free(rn);
    free(B);
    return 100.0 * nnz / np / np;
}
