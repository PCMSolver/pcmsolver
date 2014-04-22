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
 *  IntKon2.c  *
 ***************/


#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "kern.h"
#include "cubature.h"
#include "interpolate.h"
#include "intkon2.h"


#if !defined pi
#define pi 3.1415926535897932385
#endif


void IntKon2(c, element1, Q, P, M, SingleLayer, DoubleLayer, Identity)
/* GLEICHE PATCHES [0,1]^2 -> modifiziertes Skalarprodukt */
double *c;
element *element1;
cubature *Q;
vector3 ****P;
unsigned int M;
double SingleLayer(), DoubleLayer();
double Identity;
{
    unsigned int i, j;
    double w, d1, d2;
    double t1, t2, t3, t4;
    vector2 s, eta, xi, a, b;
    vector3 x, y;
    double h = 1. / (1 << element1->level);

    c[0] = c[1] = 0;
    s = vector2_make(h * element1->index_s, h * element1->index_t);
    for (i = 0; i < Q->nop; i++) {
        xi = Q->xi[i];
        w = h * h * Q->w[i] * xi.x * (1 - xi.x) * (1 - xi.x * xi.y);
        for (j = 0; j < Q->nop; j++) {
            eta = Q->xi[j];
            t1 = h * eta.x * (1 - xi.x);
            t2 = h * eta.y * (1 - xi.x * xi.y);
            t3 = t1 + h * xi.x;
            t4 = t2 + h * xi.x * xi.y;

            a.x = s.x + t1;
            a.y = s.y + t2;
            b.x = s.x + t3;
            b.y = s.y + t4;
            x = Chi(a, P[element1->patch], M);
            y = Chi(b, P[element1->patch], M);
            d1 = SingleLayer(x, y);
            d2 = DoubleLayer(x, y, n_Chi(b, P[element1->patch], M)) + DoubleLayer(y, x, n_Chi(a, P[element1->patch], M));

            a.y = s.y + t4;
            b.y = s.y + t2;
            x = Chi(a, P[element1->patch], M);
            y = Chi(b, P[element1->patch], M);
            d1 += SingleLayer(x, y);
            d2 += DoubleLayer(x, y, n_Chi(b, P[element1->patch], M)) + DoubleLayer(y, x, n_Chi(a, P[element1->patch], M));

            a.x = s.x + t2;
            a.y = s.y + t1;
            b.x = s.x + t4;
            b.y = s.y + t3;
            x = Chi(a, P[element1->patch], M);
            y = Chi(b, P[element1->patch], M);
            d1 += SingleLayer(x, y);
            d2 += DoubleLayer(x, y, n_Chi(b, P[element1->patch], M)) + DoubleLayer(y, x, n_Chi(a, P[element1->patch], M));

            a.y = s.y + t3;
            b.y = s.y + t1;
            x = Chi(a, P[element1->patch], M);
            y = Chi(b, P[element1->patch], M);
            d1 += SingleLayer(x, y);
            d2 += DoubleLayer(x, y, n_Chi(b, P[element1->patch], M)) + DoubleLayer(y, x, n_Chi(a, P[element1->patch], M));

            c[0] += 2 * w * Q->w[j] * d1;
            c[1] += w * Q->w[j] * d2;
        }
    }
    c[1] += Identity;           /* bilde +Identity */
    c[2] = c[1];
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

