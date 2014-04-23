/*
 * Purpose: Patch representation of molecular cavities from
 *			atomic coordinates and radii.
 * Version: September 21, 2009.
 * Author : Maharavo Randrianarivony.
 * Affiliation: Institut für Numerische Simulation.
 *              Rheinische Friedrich-Wilhelm Universität Bonn.
 *              Wegelerstraße 6, Bonn 53115.
 *              Germany.
 */

#include <stdlib.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "smooth.h"



void qodf_find_kehc_nerg(double **A, double *b, double *c, int n)
{
    int i, j;
    double temp;
    for (i = 0; i < n; i++) {
        temp = 0.0;
        for (j = 0; j < n; j++)
            temp = temp + A[i][j] * b[j];
        c[i] = temp;
    }
}


double sodb_lamb_fitr(double x_A, double y_A, double x_B, double y_B, double x_C, double y_C, double x, double y)
{
    int i;
    double res, **M, temp, *r, *s;
    M = (double **) malloc(3 * sizeof(double *));
    for (i = 0; i < 3; i++)
        M[i] = (double *) malloc(3 * sizeof(double));
    temp = x_A * y_B - x_A * y_C - x_B * y_A + x_B * y_C + x_C * y_A - x_C * y_B;
    M[0][0] = (y_B - y_C) / temp;
    M[0][1] = (y_C - y_A) / temp;
    M[0][2] = (y_A - y_B) / temp;
    M[1][0] = (x_C - x_B) / temp;
    M[1][1] = (x_A - x_C) / temp;
    M[1][2] = (x_B - x_A) / temp;
    M[2][0] = (x_B * y_C - x_C * y_B) / temp;
    M[2][1] = (x_C * y_A - x_A * y_C) / temp;
    M[2][2] = (x_A * y_B - x_B * y_A) / temp;
    r = (double *) malloc(3 * sizeof(double));
    s = (double *) malloc(3 * sizeof(double));
    r[0] = 1.0;
    r[1] = 0.0;
    r[2] = 0.0;
    qodf_find_kehc_nerg(M, r, s, 3);
    res = s[0] * x + s[1] * y + s[2];
    for (i = 0; i < 3; i++)
        free(M[i]);
    free(M);
    free(r);
    free(s);
    return res;
}


double cinp_lamb_rujn(double x_A, double y_A, double x_B, double y_B, double x_C, double y_C, double x, double y)
{
    int i;
    double res, **M, temp, *r, *s;
    M = (double **) malloc(3 * sizeof(double *));
    for (i = 0; i < 3; i++)
        M[i] = (double *) malloc(3 * sizeof(double));
    temp = x_A * y_B - x_A * y_C - x_B * y_A + x_B * y_C + x_C * y_A - x_C * y_B;
    M[0][0] = (y_B - y_C) / temp;
    M[0][1] = (y_C - y_A) / temp;
    M[0][2] = (y_A - y_B) / temp;
    M[1][0] = (x_C - x_B) / temp;
    M[1][1] = (x_A - x_C) / temp;
    M[1][2] = (x_B - x_A) / temp;
    M[2][0] = (x_B * y_C - x_C * y_B) / temp;
    M[2][1] = (x_C * y_A - x_A * y_C) / temp;
    M[2][2] = (x_A * y_B - x_B * y_A) / temp;
    r = (double *) malloc(3 * sizeof(double));
    s = (double *) malloc(3 * sizeof(double));
    r[0] = 0.0;
    r[1] = 1.0;
    r[2] = 0.0;
    qodf_find_kehc_nerg(M, r, s, 3);
    res = s[0] * x + s[1] * y + s[2];
    for (i = 0; i < 3; i++)
        free(M[i]);
    free(M);
    free(r);
    free(s);
    return res;
}


double goln_lamb_jocp(double x_A, double y_A, double x_B, double y_B, double x_C, double y_C, double x, double y)
{
    int i;
    double res, **M, temp, *r, *s;
    M = (double **) malloc(3 * sizeof(double *));
    for (i = 0; i < 3; i++)
        M[i] = (double *) malloc(3 * sizeof(double));
    temp = x_A * y_B - x_A * y_C - x_B * y_A + x_B * y_C + x_C * y_A - x_C * y_B;
    M[0][0] = (y_B - y_C) / temp;
    M[0][1] = (y_C - y_A) / temp;
    M[0][2] = (y_A - y_B) / temp;
    M[1][0] = (x_C - x_B) / temp;
    M[1][1] = (x_A - x_C) / temp;
    M[1][2] = (x_B - x_A) / temp;
    M[2][0] = (x_B * y_C - x_C * y_B) / temp;
    M[2][1] = (x_C * y_A - x_A * y_C) / temp;
    M[2][2] = (x_A * y_B - x_B * y_A) / temp;
    r = (double *) malloc(3 * sizeof(double));
    s = (double *) malloc(3 * sizeof(double));
    r[0] = 0.0;
    r[1] = 0.0;
    r[2] = 1.0;
    qodf_find_kehc_nerg(M, r, s, 3);
    res = s[0] * x + s[1] * y + s[2];
    for (i = 0; i < 3; i++)
        free(M[i]);
    free(M);
    free(r);
    free(s);
    return res;
}
