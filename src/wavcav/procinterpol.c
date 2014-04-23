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
#include "pln_sph.h"
#include "smooth.h"


void hiwq_find_cehq(double *t, int n, double *delta, double *sec, double *alpha, double *beta, double *gamma)
{
    int i;
    double den1, den2, a, b, *sqr;
    sqr = (double *) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        delta[i] = t[i + 1] - t[i];
        sqr[i] = delta[i] * delta[i];
    }
    for (i = 1; i < n; i++)
        sec[i] = delta[i] + delta[i - 1];
    for (i = 1; i <= n - 1; i++) {
        if (i == 1)
            den1 = sec[i];
        else
            den1 = delta[i - 2] + sec[i];
        if (i == n - 1)
            den2 = sec[i];
        else
            den2 = delta[i - 1] + sec[i + 1];

        alpha[i] = sqr[i] / den1;
        if (i == 1)
            a = delta[i] * delta[i - 1] / den1;
        else
            a = delta[i] * sec[i - 1] / den1;
        if (i == n - 1)
            b = delta[i - 1] * delta[i] / den2;
        else
            b = delta[i - 1] * sec[i + 1] / den2;
        beta[i] = a + b;
        gamma[i] = sqr[i - 1] / den2;
    }
    free(sqr);
}


void josg_lins_sung(int n, double *alpha, double *beta, double *gamma, double **MAT)
{
    int i, j;
    for (i = 0; i <= n; i++)
        for (j = 0; j <= n; j++)
            MAT[i][j] = 0.0;
    MAT[0][0] = 1.0;
    MAT[n][n] = 1.0;
    for (i = 1; i < n; i++) {
        MAT[i][i - 1] = alpha[i];
        MAT[i][i] = beta[i];
        MAT[i][i + 1] = gamma[i];
    }
}


void magr_righ_cuws(int n, double *sec, double b_f, double b_l, double *x, double *r)
{
    int i;
    r[0] = b_f;
    r[n] = b_l;
    for (i = 1; i < n; i++)
        r[i] = sec[i] * x[i];
}


void homr_find_sowk(double *t, point * X, point B_F, point B_L, int n, point * D)
{
    int i;
    double *delta, *alpha, *beta, *gamma;
    double **MAT, *rhs, *temp, *sol, *sec;

    delta = (double *) malloc(n * sizeof(double));
    alpha = (double *) malloc(n * sizeof(double));
    beta = (double *) malloc(n * sizeof(double));
    gamma = (double *) malloc(n * sizeof(double));
    sec = (double *) malloc(n * sizeof(double));
    hiwq_find_cehq(t, n, delta, sec, alpha, beta, gamma);
    MAT = (double **) malloc((n + 1) * sizeof(double *));
    for (i = 0; i < n + 1; i++)
        MAT[i] = (double *) malloc((n + 1) * sizeof(double));
    josg_lins_sung(n, alpha, beta, gamma, MAT);

    rhs = (double *) malloc((n + 1) * sizeof(double));
    temp = (double *) malloc((n + 1) * sizeof(double));
    sol = (double *) malloc((n + 1) * sizeof(double));
    for (i = 0; i <= n; i++)
        temp[i] = X[i].absi;
    magr_righ_cuws(n, sec, B_F.absi, B_L.absi, temp, rhs);
    poml_trid_fert(MAT, rhs, sol, n + 1);
    for (i = 0; i <= n; i++)
        D[i].absi = sol[i];

    for (i = 0; i <= n; i++)
        temp[i] = X[i].ordo;
    magr_righ_cuws(n, sec, B_F.ordo, B_L.ordo, temp, rhs);
    poml_trid_fert(MAT, rhs, sol, n + 1);
    for (i = 0; i <= n; i++)
        D[i].ordo = sol[i];

    for (i = 0; i <= n; i++)
        temp[i] = X[i].cote;
    magr_righ_cuws(n, sec, B_F.cote, B_L.cote, temp, rhs);
    poml_trid_fert(MAT, rhs, sol, n + 1);
    for (i = 0; i <= n; i++)
        D[i].cote = sol[i];

    free(delta);
    free(sec);
    free(temp);
    for (i = 0; i < n + 1; i++)
        free(MAT[i]);
    free(MAT);
    free(rhs);
    free(sol);
    free(alpha);
    free(beta);
    free(gamma);
}


void fumv_find_hotd(double *t, double *X, double B_F, double B_L, int n, double *D)
{
    int i;
    double *delta, *alpha, *beta, *gamma;
    double **MAT, *rhs, *temp, *sol, *sec;

    delta = (double *) malloc(n * sizeof(double));
    alpha = (double *) malloc(n * sizeof(double));
    beta = (double *) malloc(n * sizeof(double));
    gamma = (double *) malloc(n * sizeof(double));
    sec = (double *) malloc(n * sizeof(double));
    hiwq_find_cehq(t, n, delta, sec, alpha, beta, gamma);
    MAT = (double **) malloc((n + 1) * sizeof(double *));
    for (i = 0; i < n + 1; i++)
        MAT[i] = (double *) malloc((n + 1) * sizeof(double));
    josg_lins_sung(n, alpha, beta, gamma, MAT);

    rhs = (double *) malloc((n + 1) * sizeof(double));
    temp = (double *) malloc((n + 1) * sizeof(double));
    sol = (double *) malloc((n + 1) * sizeof(double));
    for (i = 0; i <= n; i++)
        temp[i] = X[i];
    magr_righ_cuws(n, sec, B_F, B_L, temp, rhs);
    poml_trid_fert(MAT, rhs, sol, n + 1);
    for (i = 0; i <= n; i++)
        D[i] = sol[i];

    free(delta);
    free(sec);
    free(temp);
    for (i = 0; i < n + 1; i++)
        free(MAT[i]);
    free(MAT);
    free(rhs);
    free(sol);
    free(alpha);
    free(beta);
    free(gamma);
}


void mofs_cubi_niws(double *t, point * X, int n, point * D)
{
    double lambda, mu;
    point B_F, B_L;
    lambda = 1.0 / 3.0;
    mu = 1.0 - lambda;
    B_F.absi = mu * X[0].absi + lambda * X[1].absi;
    B_F.ordo = mu * X[0].ordo + lambda * X[1].ordo;
    B_F.cote = mu * X[0].cote + lambda * X[1].cote;

    B_L.absi = lambda * X[n - 1].absi + mu * X[n].absi;
    B_L.ordo = lambda * X[n - 1].ordo + mu * X[n].ordo;
    B_L.cote = lambda * X[n - 1].cote + mu * X[n].cote;
    homr_find_sowk(t, X, B_F, B_L, n, D);
}


void guwn_cubi_pelc(double *t, double *X, int n, double *D)
{
    double lambda, mu;
    double B_F, B_L;
    lambda = 1.0 / 3.0;
    mu = 1.0 - lambda;
    B_F = mu * X[0] + lambda * X[1];
    B_L = lambda * X[n - 1] + mu * X[n];
    fumv_find_hotd(t, X, B_F, B_L, n, D);
}



void ralc_inte_zuts(double *t, point * X, int n, ns_curv * B)
{
    int k = 4, i;
    point *d;
    B->k = k;
    B->n = n + 2;

    for (i = 0; i < k; i++)
        B->tau[i] = t[0];
    for (i = 1; i < n; i++)
        B->tau[k + i - 1] = t[i];
    for (i = 1; i <= k; i++)
        B->tau[n + 2 + i] = t[n];

    d = (point *) malloc((n + 1) * sizeof(point));
    mofs_cubi_niws(t, X, n, d);
    for (i = 1; i <= n + 1; i++)
        getf_find_rogc_todj(d[i - 1], &B->d[i]);
    free(d);

    getf_find_rogc_todj(X[0], &B->d[0]);
    getf_find_rogc_todj(X[n], &B->d[n + 2]);

    B->prop1 = 0;
    B->prop2 = 0;
    B->prop3 = 1;
    B->prop4 = 0;
    B->v0 = t[0];
    B->v1 = t[n];
    for (i = 0; i <= n + 2; i++)
        B->w[i] = 1.0;
}


void fepn_inte_fohk(double *t, double *X, int n, ns_curv1D * B)
{
    int k = 4, i;
    double *d;
    B->k = k;
    B->n = n + 2;

    for (i = 0; i < k; i++)
        B->tau[i] = t[0];
    for (i = 1; i < n; i++)
        B->tau[k + i - 1] = t[i];
    for (i = 1; i <= k; i++)
        B->tau[n + 2 + i] = t[n];

    d = (double *) malloc((n + 1) * sizeof(double));
    guwn_cubi_pelc(t, X, n, d);
    for (i = 1; i <= n + 1; i++)
        B->d[i] = d[i - 1];
    free(d);

    B->d[0] = X[0];
    B->d[n + 2] = X[n];

    B->prop2 = 0;
    B->prop3 = 1;
    B->v0 = t[0];
    B->v1 = t[n];
    for (i = 0; i <= n + 2; i++)
        B->w[i] = 1.0;
}
