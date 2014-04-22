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

/**************
 *  sparse.c  *
 **************/


/*===================================================*
 *  Modul zur Rechnung mit duennbesetzten Matrizen.  *
 *===================================================*/


#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "sparse.h"


/*==========================================*
 *  Suchalgorithmus gemaess binary-search:  *
 *  Liefert den Index des Elements mit      *
 *  Index >= dem gesuchten Index j.         *
 *==========================================*/

unsigned int search_sparse(array, rn, j)
unsigned int *array, rn, j;
{
    unsigned int mid, low, high;

    low = 0;
    high = rn;
    while (low < high) {
        mid = (low + high) / 2;
        if (array[mid] < j)
            low = mid + 1;
        else if (array[mid] > j)
            high = mid;
        else
            return (mid);
    }
    return (low);               /* low == high! */
}


/*==================================*
 *  memmove fuer die sparse-Matrix  *
 *==================================*/

void memmove_sparse(A, i, j)
sparse *A;
unsigned int i, j;
{
    unsigned int rn, n_max;

    rn = A->row_number[i];

    if (rn == A->max_row_number[i]) {   /* Anzahl Elemente in Zeile i muss erhoeht werden */
        n_max = rn + 10;
        A->value[i] = (double *) realloc(A->value[i], n_max * sizeof(double));
        A->index[i] = (unsigned int *) realloc(A->index[i], (n_max + 1) * sizeof(unsigned int));
        A->max_row_number[i] = n_max;
    }
    memmove(&A->value[i][j + 1], &A->value[i][j], (rn - j) * sizeof(double));
    memmove(&A->index[i][j + 1], &A->index[i][j], (rn + 1 - j) * sizeof(unsigned int));
    A->row_number[i]++;
    return;
}


/*=====================================*
 *  Initialisierung der sparse-Matrix  *
 *=====================================*/

void init_sparse(A, m, n, n_max)
sparse *A;
unsigned int m, n, n_max;
{
    unsigned int i;

    A->m = m;
    A->n = n;

    A->value = (double **) malloc(m * sizeof(double *));
    A->index = (unsigned int **) malloc(m * sizeof(unsigned int *));
    A->row_number = (unsigned int *) calloc(m, sizeof(unsigned int));
    A->max_row_number = (unsigned int *) malloc(m * sizeof(unsigned int));

    for (i = 0; i < m; i++) {   /* 1 Dummy-Element: Waechter- + Pufferelement */
        A->max_row_number[i] = n_max;
        A->value[i] = (double *) malloc(n_max * sizeof(double));
        A->index[i] = (unsigned int *) malloc((n_max + 1) * sizeof(unsigned int));
        A->index[i][0] = n;     /* Waechterelement */
    }

    return;
}


/*===============================*
 *  Freigeben der sparse-Matrix  *
 *===============================*/

void free_sparse(A)
sparse *A;
{
    unsigned int i;

    for (i = 0; i < A->m; i++) {
        free(A->value[i]);
        free(A->index[i]);
    }

    free(A->value);
    free(A->index);
    free(A->row_number);
    free(A->max_row_number);

    return;
}


/*===============*
 *  A(i,j) := z  *
 *===============*/

void set_sparse(A, i, j, z)
sparse *A;
unsigned int i, j;
double z;
{
    unsigned int k;

    k = search_sparse(A->index[i], A->row_number[i], j);

    if (A->index[i][k] != j) {  /* Eintrag noch nicht vorhanden */
        memmove_sparse(A, i, k);
        A->index[i][k] = j;
    }
    A->value[i][k] = z;
    return;
}


/*===============*
 *  A(i,j) += z  *
 *===============*/

void add_sparse(A, i, j, z)
sparse *A;
unsigned int i, j;
double z;
{
    unsigned int k;

    if (z == 0)
        return;

    k = search_sparse(A->index[i], A->row_number[i], j);
    if (A->index[i][k] == j)
        A->value[i][k] += z;    /* Eintrag schon vorhanden      */
    else {                      /* Eintrag noch nicht vorhanden */
        memmove_sparse(A, i, k);
        A->index[i][k] = j;
        A->value[i][k] = z;
    }
    return;
}


/*===============*
 *  z := A(i,j)  *
 *===============*/

double get_sparse(A, i, j)
sparse *A;
unsigned int i, j;
{
    unsigned int k;

    k = search_sparse(A->index[i], A->row_number[i], j);
    if (A->index[i][k] == j)
        return (A->value[i][k]);        /* Eintrag vorhanden       */
    else
        return (0);             /* Eintrag nicht vorhanden */
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

