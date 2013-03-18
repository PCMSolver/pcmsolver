/***************
 *  IntLin1.c  *
 ***************/


/*===========================================================*
 *  Fernfeld-Quadratur-Routine:			             *
 *  Die in der Funktion init_randwerte vorab berechneten     *
 *  Auswertepunkte und Gewichte der Gauss-Quadratur werden   *	
 *  in IntLin1 zum entsprechenden Integral zusammengefuegt.  *
 *===========================================================*/


#include <stdlib.h>
#include <string.h>
#include "intvector_pwl.h"
#include "vector2.h"
#include "vector3.h"
#include "basis_pwl.h"
#include "kern.h"
#include "phi.h"
#include "cubature.h"
#include "gauss_square.h"
#include "interpolate_pwl.h"
#include "intlin1.h"
#include "integrate_pwl.h"


void IntLin1(c, element1, element2, RW, Q1, Q2, P, M, SingleLayer, DoubleLayer)
/* No-Problem-Quadrature-Routine -> modifiziertes Skalarprodukt */
double *c;
element_pwl *element1, *element2;
randwerte *RW;
cubature *Q1, *Q2;
vector3 ****P;
unsigned int M;
double SingleLayer(), DoubleLayer();
{
    double d1, d2, d3, h_s;
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
            RW->Chi[i] = Chi_pwl(xi, P[element1->patch], M);
            RW->n_Chi[i] = n_Chi_pwl(xi, P[element1->patch], M);
            RW->det_dChi[i] = h_s * Q1->w[i];
        }
    }

/* Quadratur */
    memset(c, 0, 48 * sizeof(double));
    for (i = 0; i < Q2->nop; i++) {
        eta.x = h_t * (element2->index_s + Q2->xi[i].x);
        eta.y = h_t * (element2->index_t + Q2->xi[i].y);
        y = Chi_pwl(eta, P[element2->patch], M);
        n_y = n_Chi_pwl(eta, P[element2->patch], M);
        for (j = 0; j < RW->nop; j++) {
            d1 = h_t * Q2->w[i] * RW->det_dChi[j] * SingleLayer(RW->Chi[j], y);
            d2 = h_t * Q2->w[i] * RW->det_dChi[j] * DoubleLayer(RW->Chi[j], y, n_y);
            d3 = h_t * Q2->w[i] * RW->det_dChi[j] * DoubleLayer(y, RW->Chi[j], RW->n_Chi[j]);
            Phi_times_Phi(&c[0], d1, Q1->xi[j], Q2->xi[i]);
            Phi_times_Phi(&c[16], d2, Q1->xi[j], Q2->xi[i]);
            Phi_times_Phi(&c[32], d3, Q1->xi[j], Q2->xi[i]);
        }
    }
    return;
}
