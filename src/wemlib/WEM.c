/***********
 *  WEM.c  *
 ***********/


/*==================================================*
 *  Berechnet die komprimierte Steifigkeitsmatrix.  *
 *==================================================*/


#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "sparse2.h"
#include "cubature.h"
#include "gauss_square.h"
#include "intkon1.h"
#include "integrate.h"
#include "WEM.h"


void WEM(S, W, E, T, p, M, SingleLayer, DoubleLayer, Identity)
sparse2 *S;                     /* zu berechnende Steifigkeitsmatrix            */
wavelet *W;                     /* Waveletliste                                 */
element *E;                     /* hierarchische Elementliste                   */
vector3 ****T;                  /* Koeffizienten zur Oberflaecheninterpolation  */
unsigned int p;                 /* Zahl der Patches                             */
unsigned int M;                 /* Zahl der Level                               */
double SingleLayer();           /* kernel of the single layer operator          */
double DoubleLayer();           /* kernel of the double layer operator          */
double Identity;                /* correctly scaled identity operator           */
{
    unsigned int ne;            /* Groesse von E                                */
    signed int i, j, k;         /* Laufindizes bzgl. Elemente                   */
    unsigned int *ansatz_wavelet;       /* Indizes  der Ansatz-Wavelets bzgl. Element i */
    double *weight;             /* Gewichte der Ansatz-Wavelets bzgl. Element i */
    unsigned int test_wavelet;  /* Index des Testwavelets                       */
    cubature *Q;                /* Kubatur-Formeln                              */
    randwerte *RW;              /* Randwerte                                    */
    double **prec;              /* benoetigte Quadraturgenauigkeiten            */
    double c[3];                /* berechnete Integrale Element - Element       */
    double s[3];                /* temporare Groesse fuer die Summation         */

/* Initialisierung */
    weight = NULL;
    ansatz_wavelet = NULL;
    init_randwerte(&RW, g_max);
    init_Gauss_Square(&Q, g_max);
    ne = p * (4 * (1 << 2 * M) - 1) / 3;

/* berechne log2 der Quadraturgenauigkeit */
    prec = (double **) malloc((M + 1) * sizeof(double *));
    for (i = 0; i <= M; i++) {
        prec[i] = (double *) malloc((M + 1) * sizeof(double));
        for (j = 0; j <= M; j++) {
            prec[i][j] = (2 * M - (i + j)) * (op - 2 * dp) / (2 * td + op);
            if (prec[i][j] > -fabs(i - j))
                prec[i][j] = -fabs(i - j);
            prec[i][j] += op * M - dp * (2 * M - (i + j));
        }
    }

/* berechne zeilenweise die Systemmatrix */
    for (i = 0; i < ne; i++) {
        /* Initailisierung */
        reset_randwerte(RW, g_max);
        weight = (double *) realloc(weight, E[i].wavelet_number * sizeof(double));
        ansatz_wavelet = (unsigned int *) realloc(ansatz_wavelet, E[i].wavelet_number * sizeof(unsigned int));
        memset(ansatz_wavelet, 0, E[i].wavelet_number * sizeof(unsigned int));

        /*
         * bestimme erstes Testwavelet und die Gewichte des Ansatzwavelets bzgl. Element i
         */
        test_wavelet = S->n;
        for (j = 0; j < E[i].wavelet_number; j++) {
            if (S->index[E[i].wavelet[j]][0] < test_wavelet)
                 test_wavelet = S->index[E[i].wavelet[j]][0];
            for (k = 0; W[E[i].wavelet[j]].element[k] != i; k++);
            weight[j] = W[E[i].wavelet[j]].weight[k];
        }

        /* Schleife durch die Zeilen der Systemmatrix */
        while ((test_wavelet < S->n) && (E[i].level + 2 > W[test_wavelet].level)) {
            /* Quadratur von Element i mit Wavelet test_wavelet */
            s[0] = s[1] = s[2] = 0;
            for (j = 0; j < W[test_wavelet].element_number; j++) {
                if (i > W[test_wavelet].element[j]) {
                    element_element_interaction(c, E, i, W[test_wavelet].element[j], RW, Q, T, M, prec[E[i].level][W[test_wavelet].level], SingleLayer, DoubleLayer, Identity);
                    s[0] += W[test_wavelet].weight[j] * c[0];
                    s[1] += W[test_wavelet].weight[j] * c[1];
                    s[2] += W[test_wavelet].weight[j] * c[2];
                } else if (i == W[test_wavelet].element[j]) {
                    element_element_interaction(c, E, i, W[test_wavelet].element[j], RW, Q, T, M, prec[E[i].level][W[test_wavelet].level], SingleLayer, DoubleLayer, Identity);
                    s[0] += 0.5 * W[test_wavelet].weight[j] * c[0];
                    s[1] += W[test_wavelet].weight[j] * c[1];
                }
            }

            /* Abspeichern der berechneten Integrale */
            for (j = 0; j < E[i].wavelet_number; j++) {
                if (S->index[E[i].wavelet[j]][ansatz_wavelet[j]] == test_wavelet) {
                    S->value1[E[i].wavelet[j]][ansatz_wavelet[j]] += weight[j] * s[0];
                    S->value2[E[i].wavelet[j]][ansatz_wavelet[j]] += weight[j] * s[1];
                    add_sparse2(S, test_wavelet, E[i].wavelet[j], weight[j] * s[0], weight[j] * s[2]);
                    ansatz_wavelet[j]++;
                }
            }

            /* bestimme neues Testwavelet */
            test_wavelet = S->n;
            for (j = 0; j < E[i].wavelet_number; j++) {
                if (S->index[E[i].wavelet[j]][ansatz_wavelet[j]] < test_wavelet)
                     test_wavelet = S->index[E[i].wavelet[j]][ansatz_wavelet[j]];
            }
        }

        /* loesche Element-Element-Interaktionen */
        E[i].interaction.integral_number = 0;
        free(E[i].interaction.value);
        free(E[i].interaction.index);
    }

/* Speicherplatz wieder freigeben */
    for (k = 0; k <= M; k++)
        free(prec[k]);
    free_Gauss_Square(&Q, g_max);
    free_randwerte(&RW, g_max);
    free(ansatz_wavelet);
    free(weight);
    free(prec);
}
