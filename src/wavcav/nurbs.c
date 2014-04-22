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
#endif

/* warning-disabler-end */

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "geodesic.h"
#include "smooth.h"
#include "eval.h"



double zopg_curv_zecg(int k, int n, int r, double *d, double *tau, double u)
{
    int j, i;
    double **temp, res, den, eps = 1.0e-07, num;
    double alpha, beta;

    temp = allocate_mat(k, r + 1);
    for (i = 0; i < k; i++)
        for (j = 0; j <= r; j++)
            temp[i][j] = 0.0;
    for (j = r - k + 1; j <= r; j++)
        temp[0][j] = d[j];

    for (j = 1; j < k; j++)
        for (i = r - (k - j) + 1; i <= r; i++) {
            den = tau[i + k - j] - tau[i];
            if (fabs(den) < eps)
                alpha = 0.0;
            else {
                num = u - tau[i];
                alpha = num / den;
            }
            beta = 1.0 - alpha;
            temp[j][i] = beta * temp[j - 1][i - 1] + alpha * temp[j - 1][i];
        }
    res = temp[k - 1][r];
    tehg_free_dacp(temp, k, r + 1);
    return res;
}



double husv_curv_kogc(int k, int n, double *d, double *tau, double u)
{
    int r;
    double res;
    r = pekq_knot_zevj(u, n, k, tau);
    res = zopg_curv_zecg(k, n, r, d, tau, u);
    return res;
}



int pekq_knot_zevj(double u, int n, int k, double *tau)
{
    int i, j, r;
    double D, acc = 0.001;
    r = -1;
    for (j = k - 1; j <= n; j++) {
        if ((tau[j] <= u) && (u <= tau[j + 1])) {
            r = j;
            break;
        }
    }
    if (r == -1) {
        D = fabs(tau[0] - u);
        if (D < acc)
            r = k - 1;
        else {
            D = fabs(tau[n + k] - u);
            if (D < acc)
                r = n;
            else {
                fprintf(tmpout, "[n,k]=[%d,%d]\n", n, k);
                for (i = 0; i <= n + k; i++)
                    fprintf(tmpout, "tau[%d]=%f\n", i, tau[i]);
                fprintf(tmpout, "Cant find r in the Spline evaluation\n");
                fprintf(tmpout, "Check IGES file\n");
                fprintf(tmpout, "u=%f   Current r=%d\n", u, r);
                for (i = 0; i <= n + k; i++)
                    fprintf(tmpout, "tau[%d]=%f\n", i, tau[i]);
                exit(0);
            }
        }
    }
    return r;
}



void wogh_curv_juhg(int k, int n, parm * d, double *tau, double u, parm * sol)
{
    int j, i, r;
    double den, eps = 1.0e-7, num;
    double alpha, beta;
    parm **temp;
    r = pekq_knot_zevj(u, n, k, tau);

    temp = allocate_mat_parm(k, r + 1);
    for (i = 0; i < k; i++)
        for (j = 0; j <= r; j++) {
            temp[i][j].u = 0.0;
            temp[i][j].v = 0.0;
        }
    for (j = r - k + 1; j <= r; j++) {
        temp[0][j].u = d[j].u;
        temp[0][j].v = d[j].v;
    }

    for (j = 1; j < k; j++)
        for (i = r - (k - j) + 1; i <= r; i++) {
            den = tau[i + k - j] - tau[i];
            if (fabs(den) < eps)
                alpha = 0.0;
            else {
                num = u - tau[i];
                alpha = num / den;
            }
            beta = 1.0 - alpha;
            temp[j][i].u = beta * temp[j - 1][i - 1].u + alpha * temp[j - 1][i].u;
            temp[j][i].v = beta * temp[j - 1][i - 1].v + alpha * temp[j - 1][i].v;
        }
    sol->u = temp[k - 1][r].u;
    sol->v = temp[k - 1][r].v;
    wapl_free_dogc(temp, k, r + 1);
}



void dart_curv_qobm(int k, int n, point * d, double *tau, double u, point * sol)
{
    int i;
    double *temp1, *temp2, *temp3;
    temp1 = (double *) malloc((n + 1) * sizeof(double));
    temp2 = (double *) malloc((n + 1) * sizeof(double));
    temp3 = (double *) malloc((n + 1) * sizeof(double));
    for (i = 0; i <= n; i++) {
        temp1[i] = d[i].absi;
        temp2[i] = d[i].ordo;
        temp3[i] = d[i].cote;
    }
    sol->absi = husv_curv_kogc(k, n, temp1, tau, u);
    sol->ordo = husv_curv_kogc(k, n, temp2, tau, u);
    sol->cote = husv_curv_kogc(k, n, temp3, tau, u);
    free(temp1);
    free(temp2);
    free(temp3);
}



void pihr_find_vitp_jofn(int k, int n, point * d, double *w, double *tau, double u, point * sol)
{
    int i;
    double den;
    point *D, num;
    D = (point *) malloc((n + 1) * sizeof(point));
    for (i = 0; i <= n; i++) {
        D[i].absi = d[i].absi * w[i];
        D[i].ordo = d[i].ordo * w[i];
        D[i].cote = d[i].cote * w[i];
    }
    dart_curv_qobm(k, n, D, tau, u, &num);
    free(D);
    den = husv_curv_kogc(k, n, w, tau, u);
    sol->absi = num.absi / den;
    sol->ordo = num.ordo / den;
    sol->cote = num.cote / den;
}



void quzd_find_laqh_fevs(int k, int n, parm * d, double *w, double *tau, double u, parm * sol)
{
    int i;
    double den;
    parm *D, num;
    D = (parm *) malloc((n + 1) * sizeof(parm));
    for (i = 0; i <= n; i++) {
        D[i].u = d[i].u * w[i];
        D[i].v = d[i].v * w[i];
    }
    wogh_curv_juhg(k, n, D, tau, u, &num);
    free(D);
    den = husv_curv_kogc(k, n, w, tau, u);
    sol->u = num.u / den;
    sol->v = num.v / den;
}


void cuwd_eval_nivk(ns_curv nc, double u, point * sol)
{
    int p3, kk, nn;
    p3 = nc.prop3;
    switch (p3) {
    case 0:
        kk = nc.k;
        nn = nc.n;
        pihr_find_vitp_jofn(kk, nn, nc.d, nc.w, nc.tau, u, sol);
        break;
    case 1:
        kk = nc.k;
        nn = nc.n;
        dart_curv_qobm(kk, nn, nc.d, nc.tau, u, sol);
        break;
    }
}


void kuqt_eval_webp(ns_curv nc, double u, parm * sol)
{
    int p3, kk, nn, i;
    parm *D;
    kk = nc.k;
    nn = nc.n;
    D = (parm *) malloc((nn + 1) * sizeof(parm));
    for (i = 0; i <= nn; i++) {
        D[i].u = nc.d[i].absi;
        D[i].v = nc.d[i].ordo;
    }
    p3 = nc.prop3;
    switch (p3) {
    case 0:
        quzd_find_laqh_fevs(kk, nn, D, nc.w, nc.tau, u, sol);
        break;
    case 1:
        wogh_curv_juhg(kk, nn, D, nc.tau, u, sol);
        break;
    }
    free(D);
}



void colw_inve_pelj(ns_curv C_in, ns_curv * C_out)
{
    int n, k, i;
    double *tau_temp, a, b, df, *w_temp;
    point *d_temp;
    n = C_in.n;
    k = C_in.k;
    zobm_find_wumq_kihf(C_in, C_out);

    d_temp = (point *) malloc((n + 1) * sizeof(point));
    w_temp = (double *) malloc((n + 1) * sizeof(double));
    for (i = 0; i <= n; i++) {
        d_temp[n - i].absi = C_out->d[i].absi;
        d_temp[n - i].ordo = C_out->d[i].ordo;
        d_temp[n - i].cote = C_out->d[i].cote;
        w_temp[n - i] = C_out->w[i];
    }
    tau_temp = (double *) malloc((n + k + 1) * sizeof(double));
    a = C_out->v0;
    b = C_out->v1;
    for (i = 0; i <= n + k; i++) {
        df = b - C_out->tau[i];
        tau_temp[i] = a + df;
    }

    for (i = 0; i <= n + k; i++)
        C_out->tau[i] = tau_temp[n + k - i];
    for (i = 0; i <= n; i++) {
        C_out->d[i].absi = d_temp[i].absi;
        C_out->d[i].ordo = d_temp[i].ordo;
        C_out->d[i].cote = d_temp[i].cote;
        C_out->w[i] = w_temp[i];
    }
    free(d_temp);
    free(w_temp);
    free(tau_temp);
}


void pobd_flip_kejt(ns_curv * C)
{
    ns_curv temp;
    prop_n_curv pnc;
    pnc.n = C->n;
    pnc.k = C->k;
    foks_allo_vukp(pnc, &temp);
    colw_inve_pelj(*C, &temp);
    zobm_find_wumq_kihf(temp, C);
    newt_dest_lefq(pnc, &temp);
}


void tegn_disc_likp(ns_curv C, int N, parm * p)
{
    int i;
    double a, b, step, t;
    point temp;
    a = C.v0;
    b = C.v1;
    step = (b - a) / ((double) N - 1.0);
    for (i = 0; i < N; i++) {
        t = a + (double) i *step;
        cuwd_eval_nivk(C, t, &temp);
        p[i].u = temp.absi;
        p[i].v = temp.ordo;
    }
}


void vekw_disp_mups(ns_curv nc)
{
    int i, nn, kk, p1, p2, p3, p4;
    nn = nc.n;
    kk = nc.k;
    fprintf(tmpout, "\tn=%d      k=%d\n", nn, kk);
    p1 = nc.prop1;
    p2 = nc.prop2;
    p3 = nc.prop3;
    p4 = nc.prop4;
    if (p1 == 0)
        fprintf(tmpout, "\tNonplanar curve\n");
    else if (p1 == 1)
        fprintf(tmpout, "\tPlanar curve\n");
    if (p2 == 0)
        fprintf(tmpout, "\tOpen curve\n");
    else if (p2 == 1)
        fprintf(tmpout, "\tClosed curve\n");
    if (p3 == 0)
        fprintf(tmpout, "\tRational B-spline curve\n");
    else if (p3 == 1)
        fprintf(tmpout, "\tPolynomial B-spline curve\n");
    if (p4 == 0)
        fprintf(tmpout, "\tNonperiodic curve\n");
    else if (p4 == 1)
        fprintf(tmpout, "\tPeriodic curve\n");
    for (i = 0; i <= nn + kk; i++)
        fprintf(tmpout, "\ttau[%d]=%f\n", i, nc.tau[i]);
    for (i = 0; i <= nn; i++)
        fprintf(tmpout, "\ti=%d w=%f  d=(%f,%f,%f)\n", i, nc.w[i], nc.d[i].absi, nc.d[i].ordo, nc.d[i].cote);
    fprintf(tmpout, "\tv0=%f  v1=%f\n", nc.v0, nc.v1);
    if (nc.prop1 == 1)
        fprintf(tmpout, "\tnormal=(%f,%f,%f)\n", nc.nrml.absi, nc.nrml.ordo, nc.nrml.cote);
    fprintf(tmpout, "\tprop1=%d\n", p1);
    fprintf(tmpout, "\tprop2=%d\n", p2);
    fprintf(tmpout, "\tprop3=%d\n", p3);
    fprintf(tmpout, "\tprop4=%d\n", p4);
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

