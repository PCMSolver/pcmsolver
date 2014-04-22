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

/***************
 *  IntKon1.c  *
 ***************/


/*===========================================================*
 *  Fernfeld-Quadratur-Routine:			             *
 *  Die in der Funktion init_randwerte vorab berechneten     *
 *  Auswertepunkte und Gewichte der Gauss-Quadratur werden   *
 *  in IntKon1 zum entsprechenden Integral zusammengefuegt.  *
 *===========================================================*/


#include <stdlib.h>
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "kern.h"
#include "cubature.h"
#include "gauss_square.h"
#include "interpolate.h"
#include "intkon1.h"
#include "integrate.h"


void IntKon1(c, element1, element2, RW, Q1, Q2, P, M, SingleLayer, DoubleLayer)
/* No-Problem-Quadrature-Routine -> modifiziertes Skalarprodukt */
double *c;
element *element1, *element2;
randwerte *RW;
cubature *Q1, *Q2;
vector3 ****P;
unsigned int M;
double SingleLayer(), DoubleLayer();
{
    double h_s;
    vector2 xi, eta;
    vector3 y, n_y;
    unsigned int i, j;
    double h_t = 1. / (1 << element2->level);

/* falls noetig, berechne Randwerte */
    if (RW->nop == 0) {
        RW->nop = Q1->nop;
        h_s = 1. / (1 << element1->level);
        for (i = 0; i < Q1->nop; i++) {
            xi.x = h_s * (element1->index_s + Q1->xi[i].x);
            xi.y = h_s * (element1->index_t + Q1->xi[i].y);
            RW->Chi[i] = Chi(xi, P[element1->patch], M);
            RW->n_Chi[i] = n_Chi(xi, P[element1->patch], M);
            RW->det_dChi[i] = h_s * Q1->w[i];
        }
    }
/* Quadratur */
    c[0] = c[1] = c[2] = 0;
    for (i = 0; i < Q2->nop; i++) {
        eta.x = h_t * (element2->index_s + Q2->xi[i].x);
        eta.y = h_t * (element2->index_t + Q2->xi[i].y);
        y = Chi(eta, P[element2->patch], M);
        n_y = n_Chi(eta, P[element2->patch], M);
        for (j = 0; j < RW->nop; j++) {
            c[0] += Q2->w[i] * RW->det_dChi[j] * SingleLayer(RW->Chi[j], y);
            c[1] += Q2->w[i] * RW->det_dChi[j] * DoubleLayer(RW->Chi[j], y, n_y);
            c[2] += Q2->w[i] * RW->det_dChi[j] * DoubleLayer(y, RW->Chi[j], RW->n_Chi[j]);
        }
    }

    c[0] *= h_t;                /* L^2-normiert */
    c[1] *= h_t;
    c[2] *= h_t;
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

