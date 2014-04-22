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

/*****************
 *  Integrate.c  *
 *****************/


/*=============================================================*
 *  Enthaelt alle Routinen, die im Wavelet-Galerkin-Verfahren  *
 *  fuer das aufsplitten der Integrale benoetigt werden.       *
 *=============================================================*/


#include <math.h>
#include <stdio.h>
#include "constants.h"
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "trafos.h"
#include "cubature.h"
#include "intkon1.h"
#include "intkon2.h"
#include "intkon3.h"
#include "intkon4.h"
#include "integrate.h"
#include "kern.h"


void element_element_interaction(c, E, ind1, ind2, RW, Q, R, M, prec, SingleLayer, DoubleLayer, Identity)
/* Zerlegungsalgorithmus fuer die Integration Element ind1 mit Element ind2 */
double *c;                      /* zu berechnende Integrale                      */
element *E;                     /* hierarchische Elemntliste                     */
unsigned int ind1, ind2;        /* Indizes der Integranden                       */
randwerte *RW;                  /* Randwerte bezueglich Element ind1             */
cubature *Q;                    /* Kubatur-Formeln                               */
vector3 ****R;                  /* Koeffizienten zur Oberflaecheninterpolation   */
unsigned int M;                 /* Zahl der Level                                */
double prec;                    /* Quadraturgenauigkeit                          */
double SingleLayer();           /* kernel of the single layer operator           */
double DoubleLayer();           /* kernel of the double layer operator           */
double Identity;                /* correctly scaled identity operator            */
{
    signed int j;               /* Listenindex eines Elements der Integralmatrix */
    signed int g1, g2;          /* benoetigte Quadraturgrade                     */
    double dist;                /* Abstand zwischen den Elementen ind1 und ind2  */
    unsigned int s, t;          /* Indizes der Drehungen bei sing. Integralen    */
    unsigned int CASE;          /* Fallunterscheidung bei Integration            */
    double a[12];               /* Werte der Integrale bei Unterteilung          */

/* untersuche zuerst, ob das Integral schon berechnet worden ist */
    if (ind1 >= ind2) {
        j = search_integral(&E[ind1].interaction, ind2);
        if (j != -1) {
            c[0] = E[ind1].interaction.value[j].sub[0];
            c[1] = E[ind1].interaction.value[j].sub[1];
            c[2] = E[ind1].interaction.value[j].sub[2];
            return;
        }
    } else {
        j = search_integral(&E[ind2].interaction, ind1);
        if (j != -1) {
            c[0] = E[ind2].interaction.value[j].sub[0];
            c[1] = E[ind2].interaction.value[j].sub[2];
            c[2] = E[ind2].interaction.value[j].sub[1];
            return;
        }
    }

/* bestimme den Abstand zwischen den Elementen */
    dist = distance(&E[ind1], &E[ind2]);

/* Quadratur mit Genauigkeit prec */
    if (E[ind1].level == E[ind2].level) {       /* nicht unterteilen, da beide Elemente auf gleichem Level */
        /* Falleinteilung fuer Quadratur */
        if (dist > eps)
            CASE = 1;           /* keine Gemeinsamkeiten   */
        else if (ind1 == ind2)
            CASE = 2;           /* gleiches Patch          */
        else
            CASE = compare(E[ind1].vertex, E[ind2].vertex, &s, &t);     /* gemeinsame Kante/Punkt? */

        /* Quadratur mit Genauigkeit prec */
        if (dist * (1 << E[ind1].level) < 1)
            dist = 1. / (1 << E[ind1].level);
        if (E[ind1].patch == E[ind2].patch) {   /* glatter Fall */
            quadrature_grade_smooth(&g1, &g2, E[ind1].level, E[ind2].level, dist, prec);
        } else {                /* nichtglatter Fall */
            quadrature_grade_smooth(&g1, &g2, E[ind1].level, E[ind2].level, dist, prec);
        }

        /* Wahl der Quadraturroutine aufgrund der Falleinteilung */
        switch (CASE) {
        case 1:
            IntKon1(c, &E[ind1], &E[ind2], &RW[g1], &Q[g1], &Q[g1], R, M, SingleLayer, DoubleLayer);
            break;
        case 2:
            IntKon2(c, &E[ind1], &Q[g1], R, M, SingleLayer, DoubleLayer, Identity);
            break;
        case 3:
            IntKon3(c, &E[ind1], &E[ind2], s, t, &Q[g1], R, M, SingleLayer, DoubleLayer);
            break;
        case 4:
            IntKon4(c, &E[ind1], &E[ind2], s, t, &Q[g1], R, M, SingleLayer, DoubleLayer);
            break;
        default:
            c[0] = c[1] = c[2] = 0;
            printf("ERROR: CASE != 1,2,3,4 \n");
        }
    } else {                    /* (E[ind1].level > E[ind2].level) */
        if (dist * (1 << E[ind2].level) >= q) { /* Quadratur */
            quadrature_grade_smooth(&g1, &g2, E[ind1].level, E[ind2].level, dist, prec);
            IntKon1(c, &E[ind1], &E[ind2], &RW[g1], &Q[g1], &Q[g2], R, M, SingleLayer, DoubleLayer);
        } else {                /* adaptives Unterteilen */
            element_element_interaction(&a[0], E, ind1, E[ind2].son[0], RW, Q, R, M, prec, SingleLayer, DoubleLayer, Identity);
            element_element_interaction(&a[3], E, ind1, E[ind2].son[1], RW, Q, R, M, prec, SingleLayer, DoubleLayer, Identity);
            element_element_interaction(&a[6], E, ind1, E[ind2].son[2], RW, Q, R, M, prec, SingleLayer, DoubleLayer, Identity);
            element_element_interaction(&a[9], E, ind1, E[ind2].son[3], RW, Q, R, M, prec, SingleLayer, DoubleLayer, Identity);

            c[0] = 0.5 * (a[0] + a[3] + a[6] + a[9]);
            c[1] = 0.5 * (a[1] + a[4] + a[7] + a[10]);
            c[2] = 0.5 * (a[2] + a[5] + a[8] + a[11]);
        }
    }

/* Abspeichern des Integrals */
    if (ind1 >= ind2)
        set_integral(&E[ind1].interaction, ind2, c);
    else {
        a[0] = c[0];
        a[1] = c[2];
        a[2] = c[1];
        set_integral(&E[ind2].interaction, ind1, a);
    }
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

