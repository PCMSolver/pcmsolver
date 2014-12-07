/***************
 *  sparse2.c  *
 ***************/


/*===================================================*
 *  Modul zur Rechnung mit duennbesetzten Matrizen.  *
 *===================================================*/


#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "sparse2.h"


/*==========================================*
 *  Suchalgorithmus gemaess binary-search:  *
 *  Liefert den Index des Elements mit      *
 *  Index >= dem gesuchten Index j.         *
 *==========================================*/

unsigned int search_sparse2(array, rn, j)
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

void memmove_sparse2(A, i, j)
sparse2 *A;
unsigned int i, j;
{
    unsigned int rn, n_max;

    rn = A->row_number[i];

    if (rn == A->max_row_number[i]) {   /* Anzahl Elemente in Zeile i muss erhoeht werden */
        n_max = rn + 10;
        A->value1[i] = (double *) realloc(A->value1[i], n_max * sizeof(double));
        A->value2[i] = (double *) realloc(A->value2[i], n_max * sizeof(double));
        A->index[i] = (unsigned int *) realloc(A->index[i], (n_max + 1) * sizeof(unsigned int));
        A->max_row_number[i] = n_max;
    }
    memmove(&A->value1[i][j + 1], &A->value1[i][j], (rn - j) * sizeof(double));
    memmove(&A->value2[i][j + 1], &A->value2[i][j], (rn - j) * sizeof(double));
    memmove(&A->index[i][j + 1], &A->index[i][j], (rn + 1 - j) * sizeof(unsigned int));
    A->row_number[i]++;
    return;
}


/*===================================*
 *  memmove fuer die Matrix-Pattern  *
 *===================================*/

void memmove_pattern(A, i, j)
sparse2 *A;
unsigned int i, j;
{
    unsigned int rn, n_max;

    rn = A->row_number[i];

    if (rn == A->max_row_number[i]) {   /* Anzahl Elemente in Zeile i muss erhoeht werden */
        n_max = rn + 10;
        A->index[i] = (unsigned int *) realloc(A->index[i], (n_max + 1) * sizeof(unsigned int));
        A->max_row_number[i] = n_max;
    }
    memmove(&A->index[i][j + 1], &A->index[i][j], (rn + 1 - j) * sizeof(unsigned int));
    A->row_number[i]++;
    return;
}


/*=====================================*
 *  Initialisierung der sparse-Matrix  *
 *=====================================*/

void init_sparse2(A, m, n, n_max)
sparse2 *A;
unsigned int m, n, n_max;
{
    unsigned int i;

    A->m = m;
    A->n = n;

    A->value1 = (double **) malloc(m * sizeof(double *));
    A->value2 = (double **) malloc(m * sizeof(double *));
    A->index = (unsigned int **) malloc(m * sizeof(unsigned int *));
    A->row_number = (unsigned int *) calloc(m, sizeof(unsigned int));
    A->max_row_number = (unsigned int *) malloc(m * sizeof(unsigned int));

    for (i = 0; i < m; i++) {   /* 1 Dummy-Element: Waechter- + Pufferelement */
        A->max_row_number[i] = n_max;
        A->value1[i] = (double *) malloc(n_max * sizeof(double));
        A->value2[i] = (double *) malloc(n_max * sizeof(double));
        A->index[i] = (unsigned int *) malloc((n_max + 1) * sizeof(unsigned int));
        A->index[i][0] = n;     /* Waechterelement */
    }

    return;
}


/*======================================*
 *  Initialisierung der Matrix-Pattern  *
 *======================================*/

void init_pattern(A, m, n, n_max)
sparse2 *A;
unsigned int m, n, n_max;
{
    unsigned int i;

    A->m = m;
    A->n = n;

    A->value1 = (double **) malloc(m * sizeof(double *));
    A->value2 = (double **) malloc(m * sizeof(double *));
    A->index = (unsigned int **) malloc(m * sizeof(unsigned int *));
    A->row_number = (unsigned int *) calloc(m, sizeof(unsigned int));
    A->max_row_number = (unsigned int *) malloc(m * sizeof(unsigned int));

    for (i = 0; i < m; i++) {   /* 1 Dummy-Element: Waechter- + Pufferelement */
        A->max_row_number[i] = n_max;
        A->index[i] = (unsigned int *) malloc((n_max + 1) * sizeof(unsigned int));
        A->index[i][0] = n;     /* Waechterelement */
    }

    return;
}


/*===============================*
 *  Freigeben der sparse-Matrix  *
 *===============================*/

void free_sparse2(A)
sparse2 *A;
{
    unsigned int i;

    for (i = 0; i < A->m; i++) {
        free(A->index[i]);
        free(A->value1[i]);
        free(A->value2[i]);
    }

    free(A->index);
    free(A->value1);
    free(A->value2);
    free(A->row_number);
    free(A->max_row_number);

    return;
}


/*==================================*
 *  Fuege A(i,j) den Pattern hinzu  *
 *==================================*/

void set_pattern(A, i, j)
sparse2 *A;
unsigned int i, j;
{
    unsigned int k;

    k = search_sparse2(A->index[i], A->row_number[i], j);

    if (A->index[i][k] != j) {  /* Eintrag noch nicht vorhanden */
        memmove_pattern(A, i, k);
        A->index[i][k] = j;
    }
    return;
}


/*=====================*
 *  A(i,j) := (z1,z2)  *
 *=====================*/

void set_sparse2(A, i, j, z1, z2)
sparse2 *A;
unsigned int i, j;
double z1, z2;
{
    unsigned int k;

    k = search_sparse2(A->index[i], A->row_number[i], j);

    if (A->index[i][k] != j) {  /* Eintrag noch nicht vorhanden */
        memmove_sparse2(A, i, k);
        A->index[i][k] = j;
    }
    A->value1[i][k] = z1;
    A->value2[i][k] = z2;
    return;
}


/*=====================*
 *  A(i,j) += (z1,z2)  *
 *=====================*/

void add_sparse2(A, i, j, z1, z2)
sparse2 *A;
unsigned int i, j;
double z1, z2;
{
    unsigned int k;

    if ((z1 == 0) && (z2 == 0))
        return;

    k = search_sparse2(A->index[i], A->row_number[i], j);
    if (A->index[i][k] == j) {  /* Eintrag schon vorhanden      */
        A->value1[i][k] += z1;
        A->value2[i][k] += z2;
    } else {                    /* Eintrag noch nicht vorhanden */
        memmove_sparse2(A, i, k);
        A->value1[i][k] = z1;
        A->value2[i][k] = z2;
        A->index[i][k] = j;
    }
    return;
}


/*=====================*
 *  (z1,z2) := A(i,j)  *
 *=====================*/

void get_sparse2(A, i, j, z1, z2)
sparse2 *A;
unsigned int i, j;
double *z1, *z2;
{
    unsigned int k;

    k = search_sparse2(A->index[i], A->row_number[i], j);
    if (A->index[i][k] == j) {  /* Eintrag vorhanden       */
        *z1 = A->value1[i][k];
        *z2 = A->value2[i][k];
    } else
        *z1 = *z2 = 0;          /* Eintrag nicht vorhanden */
    return;
}

/*====================*
 *  nz = non_zero(A)  *
 *====================*/

unsigned int non_zero(A)
sparse2 *A;
{
    unsigned int i, nz;

    nz = 0;
    for (i = 0; i < A->m; i++)
        nz += A->row_number[i];

    return (nz);
}


/*=====================*
 *  Garbage-Collector  *
 *=====================*/

void garbage_collect(A)
sparse2 *A;
{
    unsigned int i, rn;

    for (i = 0; i < A->m; i++) {
        rn = A->row_number[i];
        A->value1[i] = (double *) realloc(A->value1[i], rn * sizeof(double));
        A->value2[i] = (double *) realloc(A->value2[i], rn * sizeof(double));
        A->index[i] = (unsigned int *) realloc(A->index[i], (rn + 1) * sizeof(unsigned int));
        A->max_row_number[i] = rn;
    }

    return;
}


/*==================*
 *  Finish-Pattern  *
 *==================*/

void finish_pattern(A)
sparse2 *A;
{
    unsigned int i, rn;

    for (i = 0; i < A->m; i++) {
        rn = A->row_number[i];
        A->value1[i] = (double *) calloc(rn, sizeof(double));
        A->value2[i] = (double *) calloc(rn, sizeof(double));
        A->index[i] = (unsigned int *) realloc(A->index[i], (rn + 1) * sizeof(unsigned int));
        A->max_row_number[i] = rn;
    }

    return;
}

void print_sparse2(sparse2 *A)
{
    int i, j;
    for (i = 0; i < A->m; i++) {
        for (j = 0; j < A->max_row_number[i]; j++) {
            printf("%d %d %d %18.6f %18.6f\n", i, j, A->index[i][j], A->value1[i][j], A->value2[i][j]);
        }
    }
}

void fprint_sparse2(sparse2 *A, char * infile)
{
    FILE *fp = fopen(infile, "w");
    int i, j;
    for (i = 0; i < A->m; i++) {
        for (j = 0; j < A->max_row_number[i]; j++) {
            fprintf(fp, "%d %d %d %18.6f %18.6f\n", i, j, A->index[i][j], A->value1[i][j], A->value2[i][j]);
        }
    }
    fclose(fp);
}

void fprint_vec(double *A, int elements, char * infile)
{
    FILE *fp = fopen(infile, "w");
    for (int i = 0; i < elements; i++) {
        fprintf(fp, "%d %18.6f\n", i, A[i]);
    }
    fclose(fp);
}

