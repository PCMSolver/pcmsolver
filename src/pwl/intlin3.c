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
 *  IntLin3.c  *
 ***************/


#include <string.h>
#include "intvector_pwl.h"
#include "vector2.h"
#include "vector3.h"
#include "basis_pwl.h"
#include "trafos_pwl.h"
#include "kern.h"
#include "phi.h"
#include "interpolate_pwl.h"
#include "cubature.h"
#include "intlin3.h"


void IntLin3(c, element1, element2, ind_s, ind_t, Q, P, M, SingleLayer, DoubleLayer)
/* GEMEINSAME KANTE [0,1] -> modifiziertes Skalarprodukt */
double *c;
element_pwl *element1, *element2;
unsigned int ind_s, ind_t;
cubature *Q;
vector3 ****P;
unsigned int M;
double SingleLayer(), DoubleLayer();
{
    unsigned int i, j;
    double d1, d2, d3, w, t1, t2, t3, t4;
    vector2 s, t, xi, eta, a, b, u, v;
    vector3 x, y;
    double h = 1. / (1 << element1->level);

    memset(c, 0, 48 * sizeof(double));
    s = vector2_make(h * element1->index_s, h * element1->index_t);
    t = vector2_make(h * element2->index_s, h * element2->index_t);

    for (i = 0; i < Q->nop; i++) {
        xi = Q->xi[i];
        w = h * h * xi.y * xi.y * Q->w[i];
        t1 = xi.x * (1 - xi.y);
        t2 = (1 - xi.x) * (1 - xi.y);

        for (j = 0; j < Q->nop; j++) {
            eta = vector2_Smul(xi.y, Q->xi[j]);
            t3 = xi.x * (1 - eta.x);
            t4 = (1 - xi.x) * (1 - eta.x);

            a = Tau_pwl(t1, eta.x, ind_s);
            b = Tau_pwl(t2, eta.y, ind_t);
            u = Kappa_pwl(s, a, h);
            v = Kappa_pwl(t, b, h);
            x = Chi_pwl(u, P[element1->patch], M);
            y = Chi_pwl(v, P[element2->patch], M);
            d1 = w * Q->w[j] * (1 - xi.y) * SingleLayer(x, y);
            d2 = w * Q->w[j] * (1 - xi.y) * DoubleLayer(x, y, n_Chi_pwl(v, P[element2->patch], M));
            d3 = w * Q->w[j] * (1 - xi.y) * DoubleLayer(y, x, n_Chi_pwl(u, P[element1->patch], M));
            Phi_times_Phi(&c[0], d1, a, b);
            Phi_times_Phi(&c[16], d2, a, b);
            Phi_times_Phi(&c[32], d3, a, b);

            a = Tau_pwl(1 - t1, eta.x, ind_s);
            b = Tau_pwl(1 - t2, eta.y, ind_t);
            u = Kappa_pwl(s, a, h);
            v = Kappa_pwl(t, b, h);
            x = Chi_pwl(u, P[element1->patch], M);
            y = Chi_pwl(v, P[element2->patch], M);
            d1 = w * Q->w[j] * (1 - xi.y) * SingleLayer(x, y);
            d2 = w * Q->w[j] * (1 - xi.y) * DoubleLayer(x, y, n_Chi_pwl(v, P[element2->patch], M));
            d3 = w * Q->w[j] * (1 - xi.y) * DoubleLayer(y, x, n_Chi_pwl(u, P[element1->patch], M));
            Phi_times_Phi(&c[0], d1, a, b);
            Phi_times_Phi(&c[16], d2, a, b);
            Phi_times_Phi(&c[32], d3, a, b);

            a = Tau_pwl(t3, xi.y, ind_s);
            b = Tau_pwl(t4, eta.y, ind_t);
            u = Kappa_pwl(s, a, h);
            v = Kappa_pwl(t, b, h);
            x = Chi_pwl(u, P[element1->patch], M);
            y = Chi_pwl(v, P[element2->patch], M);
            d1 = w * Q->w[j] * (1 - eta.x) * SingleLayer(x, y);
            d2 = w * Q->w[j] * (1 - eta.x) * DoubleLayer(x, y, n_Chi_pwl(v, P[element2->patch], M));
            d3 = w * Q->w[j] * (1 - eta.x) * DoubleLayer(y, x, n_Chi_pwl(u, P[element1->patch], M));
            Phi_times_Phi(&c[0], d1, a, b);
            Phi_times_Phi(&c[16], d2, a, b);
            Phi_times_Phi(&c[32], d3, a, b);

            a = Tau_pwl(1 - t3, xi.y, ind_s);
            b = Tau_pwl(1 - t4, eta.y, ind_t);
            u = Kappa_pwl(s, a, h);
            v = Kappa_pwl(t, b, h);
            x = Chi_pwl(u, P[element1->patch], M);
            y = Chi_pwl(v, P[element2->patch], M);
            d1 = w * Q->w[j] * (1 - eta.x) * SingleLayer(x, y);
            d2 = w * Q->w[j] * (1 - eta.x) * DoubleLayer(x, y, n_Chi_pwl(v, P[element2->patch], M));
            d3 = w * Q->w[j] * (1 - eta.x) * DoubleLayer(y, x, n_Chi_pwl(u, P[element1->patch], M));
            Phi_times_Phi(&c[0], d1, a, b);
            Phi_times_Phi(&c[16], d2, a, b);
            Phi_times_Phi(&c[32], d3, a, b);

            a = Tau_pwl(t4, eta.y, ind_s);
            b = Tau_pwl(t3, xi.y, ind_t);
            u = Kappa_pwl(s, a, h);
            v = Kappa_pwl(t, b, h);
            x = Chi_pwl(u, P[element1->patch], M);
            y = Chi_pwl(v, P[element2->patch], M);
            d1 = w * Q->w[j] * (1 - eta.x) * SingleLayer(x, y);
            d2 = w * Q->w[j] * (1 - eta.x) * DoubleLayer(x, y, n_Chi_pwl(v, P[element2->patch], M));
            d3 = w * Q->w[j] * (1 - eta.x) * DoubleLayer(y, x, n_Chi_pwl(u, P[element1->patch], M));
            Phi_times_Phi(&c[0], d1, a, b);
            Phi_times_Phi(&c[16], d2, a, b);
            Phi_times_Phi(&c[32], d3, a, b);

            a = Tau_pwl(1 - t4, eta.y, ind_s);
            b = Tau_pwl(1 - t3, xi.y, ind_t);
            u = Kappa_pwl(s, a, h);
            v = Kappa_pwl(t, b, h);
            x = Chi_pwl(u, P[element1->patch], M);
            y = Chi_pwl(v, P[element2->patch], M);
            d1 = w * Q->w[j] * (1 - eta.x) * SingleLayer(x, y);
            d2 = w * Q->w[j] * (1 - eta.x) * DoubleLayer(x, y, n_Chi_pwl(v, P[element2->patch], M));
            d3 = w * Q->w[j] * (1 - eta.x) * DoubleLayer(y, x, n_Chi_pwl(u, P[element1->patch], M));
            Phi_times_Phi(&c[0], d1, a, b);
            Phi_times_Phi(&c[16], d2, a, b);
            Phi_times_Phi(&c[32], d3, a, b);
        }
    }
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

