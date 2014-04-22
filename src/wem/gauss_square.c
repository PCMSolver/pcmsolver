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

/********************
 *  Gauss_Square.c  *
 ********************/


/*============================================*
 *  Definiert die Gauss-Quadraturformeln auf  *
 *  dem Referenzviereck [0,1]^2.              *
 *============================================*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vector2.h"
#include "quadrature.h"
#include "gauss_legendre.h"
#include "cubature.h"
#include "gauss_square.h"


void Tensor_Regel_Square(Q, R)
/* bastelt die Tensor-Produkt-Regel zusammen */
cubature *Q;
quadrature *R;
{
    unsigned int k;

    Q->nop = R->nop * R->nop;
    Q->xi = (vector2 *) malloc(Q->nop * sizeof(vector2));
    Q->w = (double *) malloc(Q->nop * sizeof(double));

    for (k = 0; k < Q->nop; k++) {
        Q->xi[k] = vector2_make(R->xi[k / R->nop], R->xi[k % R->nop]);
        Q->w[k] = R->w[k / R->nop] * R->w[k % R->nop];
    }
    return;
}


void init_Gauss_Square(Q, g)
cubature **Q;
unsigned int g;
{
    quadrature *R;
    unsigned int k;

/* Fehler-Routine */
    if (g > 15) {
        printf("g should be less equal 15\n");
        exit(0);
    }
/* Tensor-Produkte zusammenbasteln */
    init_Gauss_Legendre(&R, g);
    (*Q) = (cubature *) malloc(g * sizeof(cubature));
    for (k = 0; k < g; k++)
        Tensor_Regel_Square(&(*Q)[k], &R[k]);
    free_Gauss_Legendre(&R, g);

    return;
}


void free_Gauss_Square(Q, g)
cubature **Q;
unsigned int g;
{
    unsigned int k;

    for (k = 0; k < g; k++) {
        free((*Q)[k].xi);
        free((*Q)[k].w);
    }
    free(*Q);

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

