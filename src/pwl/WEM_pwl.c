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
 *  WEM_pwl.c  *
 ***********/


/*==================================================*
 *  Berechnet die komprimierte Steifigkeitsmatrix.  *
 *==================================================*/


#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "intvector_pwl.h"
#include "vector2.h"
#include "vector3.h"
#include "basis_pwl.h"
#include "sparse2.h"
#include "cubature.h"
#include "gauss_square.h"
#include "intlin1.h"
#include "integrate_pwl.h"
#include "WEM_pwl.h"


void WEM_pwl(S, W, P, E, T, p, M, SingleLayer, DoubleLayer, Identity)
sparse2 *S;                     /* zu berechnende Steifigkeitsmatrix            */
wavelet_pwl *W;                 /* Waveletliste                                 */
vector3 *P;                     /* Punkteliste                                  */
element_pwl *E;                 /* hierarchische Elementliste                   */
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
    unsigned int *weight_index; /* Gewichte der Ansatz-Wavelets bzgl. Element i */
    unsigned int test_wavelet;  /* Index des Testwavelets                       */
    cubature *Q;                /* Kubatur-Formeln                              */
    randwerte *RW;              /* Randwerte                                    */
    double **prec;              /* benoetigte Quadraturgenauigkeiten            */
    double c[48];               /* berechnete Integrale Element - Element       */
    double s[12], t[3];         /* temporare Groessen fuer die Summation        */

/* Initialisierung */
    weight_index = NULL;
    ansatz_wavelet = NULL;
    init_randwerte(&RW, g_max);
    init_Gauss_Square(&Q, g_max);
    ne = p * (4 * (1 << 2 * M) - 1) / 3;

/* berechne log2 der Quadraturgenauigkeit */
    prec = (double **) malloc((M + 1) * sizeof(double *));
    for (i = 0; i <= M; i++) {
        prec[i] = (double *) malloc((M + 1) * sizeof(double));
        for (j = 0; j <= M; j++) {
            prec[i][j] = (2 * M - (i + j)) * (op - 2 * dp) / (2 * td_pwl + op);
            if (prec[i][j] > -fabs(i - j))
                prec[i][j] = -fabs(i - j);
            prec[i][j] += op * M - dp * (2 * M - (i + j));
        }
    }

/* berechne zeilenweise die Systemmatrix */
    for (i = 0; i < ne; i++) {
        /* Initailisierung */
        reset_randwerte(RW, g_max);
        weight_index = (unsigned int *) realloc(weight_index, E[i].wavelet_number * sizeof(unsigned int));
        ansatz_wavelet = (unsigned int *) realloc(ansatz_wavelet, E[i].wavelet_number * sizeof(unsigned int));
        memset(ansatz_wavelet, 0, E[i].wavelet_number * sizeof(unsigned int));

        /* bestimme erstes Testwavelet und die Gewichte des Ansatzwavelets bzgl. Element i */
        test_wavelet = S->n;
        for (j = 0; j < E[i].wavelet_number; j++) {
            if (S->index[E[i].wavelet[j]][0] < test_wavelet)
                 test_wavelet = S->index[E[i].wavelet[j]][0];
            for (weight_index[j] = 0; W[E[i].wavelet[j]].element[weight_index[j]] != i; weight_index[j]++);
        }

        /* Schleife durch die Zeilen der Systemmatrix */
        while ((test_wavelet < S->n) && (E[i].level + 2 > W[test_wavelet].level)) {
            /* Quadratur von Element i mit Wavelet test_wavelet */
            memset(s, 0, 12 * sizeof(double));
            for (j = 0; j < W[test_wavelet].element_number; j++) {
                if (i > W[test_wavelet].element[j]) {
                    element_element_interaction_pwl(c, P, E, i, W[test_wavelet].element[j], RW, Q, T, M, 
                                                    prec[E[i].level][W[test_wavelet].level], SingleLayer, 
                                                    DoubleLayer, Identity);
                    for (k = 0; k < 48; k++)
                        s[k / 4] += W[test_wavelet].weight[j][k % 4] * c[k];
                } else if (i == W[test_wavelet].element[j]) {
                    element_element_interaction_pwl(c, P, E, i, W[test_wavelet].element[j], RW, Q, T, M, 
                                                    prec[E[i].level][W[test_wavelet].level], SingleLayer, 
                                                    DoubleLayer, Identity);
                    for (k = 0; k < 48; k++)
                        s[k / 4] += 0.5 * W[test_wavelet].weight[j][k % 4] * c[k];
                }
            }

            /* Abspeichern der berechneten Integrale */
            for (j = 0; j < E[i].wavelet_number; j++) {
                if (S->index[E[i].wavelet[j]][ansatz_wavelet[j]] == test_wavelet) {
                    t[0] = t[1] = t[2] = 0;
                    for (k = 0; k < 12; k++)
                        t[k / 4] += W[E[i].wavelet[j]].weight[weight_index[j]][k % 4] * s[k];
                    S->value1[E[i].wavelet[j]][ansatz_wavelet[j]] += t[0];
                    S->value2[E[i].wavelet[j]][ansatz_wavelet[j]] += t[1];
                    add_sparse2(S, test_wavelet, E[i].wavelet[j], t[0], t[2]);
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
    free(weight_index);
    free(prec);
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

