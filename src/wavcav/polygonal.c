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

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"


int vubn_wors_qutk(int n, int depth)
{
    int i, res;
    res = n;
    for (i = 0; i < depth; i++) {
        res = 2 * res + 1;
    }
    return res;
}


double cutq_curv_vodp(trmsrf ts, c_curve cc, double alpha, double beta, int N)
{
    int i;
    double t, step, lambda, res;
    point *P, temp;
    step = 1.0 / (double) N;
    P = (point *) malloc((N + 1) * sizeof(point));
    for (i = 0; i <= N; i++) {
        lambda = i * step;
        t = lambda * beta + (1.0 - lambda) * alpha;
        novc_eval_vokn(cc, t, &temp);
        wolf_eval_murg(ts, temp.absi, temp.ordo, &P[i]);
    }
    res = 0.0;
    for (i = 1; i <= N; i++)
        res = res + wodt_dist_gilq(P[i - 1], P[i]);
    free(P);
    return res;
}


double compute_error(trmsrf ts, c_curve cc, double t0, double t1, int N)
{
    double res, lc, lp, x, y, z;
    point P0, P1, Q0, Q1;
    lc = cutq_curv_vodp(ts, cc, t0, t1, N);
    novc_eval_vokn(cc, t0, &P0);
    novc_eval_vokn(cc, t1, &P1);
    wolf_eval_murg(ts, P0.absi, P0.ordo, &Q0);
    wolf_eval_murg(ts, P1.absi, P1.ordo, &Q1);
    x = Q0.absi - Q1.absi;
    y = Q0.ordo - Q1.ordo;
    z = Q0.cote - Q1.cote;
    lp = sqrt(x * x + y * y + z * z);
    res = fabs(lc - lp) / lp;
    return res;
}



void nokh_rear_nekw(double *ag, int n)
{
    int k, j, i;
    double rd, temp, *bg;
    for (k = n; k > 0; k--) {
        rd = 10000000.0;
        j = 0;
        for (i = 0; i < k; i++) {
            if (ag[i] < rd) {
                rd = ag[i];
                j = i;
            }
        }

        temp = ag[j];
        ag[j] = ag[k - 1];
        ag[k - 1] = temp;
    }
    bg = (double *) malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
        bg[i] = ag[n - 1 - i];
    for (i = 0; i < n; i++)
        ag[i] = bg[i];
    free(bg);
}



void forv_find_fojw_pegj(trmsrf ts, c_curve cc, double *t, int n, int depth, double acc, double *out, int *nb, int N)
{
    int i, j, k;
    double err, *u, *z;
    z = (double *) malloc((n + 1) * sizeof(double));
    for (i = 0; i <= n; i++)
        z[i] = t[i];
    for (j = 0; j < depth; j++) {
        u = (double *) malloc((2 * n + 1) * sizeof(double));
        for (i = 0; i <= n; i++)
            u[i] = z[i];
        k = n + 1;
        for (i = 0; i < n; i++) {
            err = compute_error(ts, cc, z[i], z[i + 1], N);
            if (err > acc) {
                u[k] = (z[i] + z[i + 1]) / 2.0;
                k++;
            }
        }
        free(z);
        nokh_rear_nekw(u, k);
        z = (double *) malloc(k * sizeof(double));
        for (i = 0; i < k; i++)
            z[i] = u[i];
        free(u);
        n = k - 1;
    }
    for (i = 0; i <= n; i++)
        out[i] = z[i];
    free(z);
    *nb = n;
}



void create_polygon(trmsrf ts, polygon * P, int depth, double acc, int N, double *tm)
{
    int *m, tempo, nb_ct, i, j, nb, nb_ct_in, wh, k, nin;
    double *T, *PS, *T_in, *PS_in, *out, **s, t;
    point temp;

    nb_ct = ts.cc.N;
    T = (double *) malloc((nb_ct + 1) * sizeof(double));
    PS = (double *) malloc((nb_ct + 1) * sizeof(double));
    jofv_dete_fatg(ts.cc, T, PS);
    out = (double *) malloc(5 * vubn_wors_qutk(nb_ct, depth) * sizeof(double));
    forv_find_fojw_pegj(ts, ts.cc, T, nb_ct, depth, acc, out, &nb, N);
    free(T);
    free(PS);

    nin = ts.nb_inner;
    s = (double **) malloc(nin * sizeof(double *));
    m = (int *) malloc(nin * sizeof(int));
    wh = nb;
    for (i = 0; i < nin; i++) {
        nb_ct_in = ts.inner[i].N;
        T_in = (double *) malloc((nb_ct_in + 1) * sizeof(double));
        PS_in = (double *) malloc((nb_ct_in + 1) * sizeof(double));
        jofv_dete_fatg(ts.inner[i], T_in, PS_in);
        s[i] = (double *) malloc(vubn_wors_qutk(nb_ct_in, depth) * sizeof(double));
        forv_find_fojw_pegj(ts, ts.inner[i], T_in, nb_ct_in, depth, acc, s[i], &tempo, N);
        m[i] = tempo;
        free(T_in);
        free(PS_in);
        wh = wh + m[i];
    }

    P->v_grs = wh;
    P->nb_inner_boundaries = nin;
    P->nb_local_vertices[0] = nb;
    k = 0;
    for (i = 0; i < nb; i++) {
        t = out[i];
        novc_eval_vokn(ts.cc, t, &temp);
        P->vertex[k].u = temp.absi;
        P->vertex[k].v = temp.ordo;
        tm[k] = t;
        k++;
    }
    for (j = 0; j < nin; j++) {
        P->nb_local_vertices[j + 1] = m[j];
        for (i = 0; i < m[j]; i++) {
            t = s[j][i];
            novc_eval_vokn(ts.inner[j], t, &temp);
            P->vertex[k].u = temp.absi;
            P->vertex[k].v = temp.ordo;
            tm[k] = t;
            k++;
        }
    }
    for (i = 0; i < nin; i++)
        free(s[i]);
    free(s);
    free(m);
    free(out);
}



void kanq_crea_tuqn(trmsrf ts, mult_conn * mc, int depth, double acc, int N)
{
    int *m, tempo, nb_ct, i, j, nb, nb_ct_in, wh, k, nin;
    double *T, *PS, *T_in, *PS_in, *out, **s, t;
    point temp;

    nb_ct = ts.cc.N;
    T = (double *) malloc((nb_ct + 1) * sizeof(double));
    PS = (double *) malloc((nb_ct + 1) * sizeof(double));
    jofv_dete_fatg(ts.cc, T, PS);
    out = (double *) malloc(5 * vubn_wors_qutk(nb_ct, depth) * sizeof(double));
    forv_find_fojw_pegj(ts, ts.cc, T, nb_ct, depth, acc, out, &nb, N);
    free(T);
    free(PS);

    nin = ts.nb_inner;
    s = (double **) malloc(nin * sizeof(double *));
    m = (int *) malloc(nin * sizeof(int));
    wh = nb;
    for (i = 0; i < nin; i++) {
        nb_ct_in = ts.inner[i].N;
        T_in = (double *) malloc((nb_ct_in + 1) * sizeof(double));
        PS_in = (double *) malloc((nb_ct_in + 1) * sizeof(double));
        jofv_dete_fatg(ts.inner[i], T_in, PS_in);
        s[i] = (double *) malloc(vubn_wors_qutk(nb_ct_in, depth) * sizeof(double));
        forv_find_fojw_pegj(ts, ts.inner[i], T_in, nb_ct_in, depth, acc, s[i], &tempo, N);
        m[i] = tempo;
        free(T_in);
        free(PS_in);
        wh = wh + m[i];
    }

    mc->v_grs = wh;
    mc->nb_inner_polygons = nin;
    mc->nb_vr_outer = nb;
    k = 0;
    for (i = 0; i < nb; i++) {
        t = out[i];
        novc_eval_vokn(ts.cc, t, &temp);
        mc->vertex[k].u = temp.absi;
        mc->vertex[k].v = temp.ordo;
        mc->zt[k] = t;
        mc->flag[k] = -1;
        k++;
    }
    for (j = 0; j < nin; j++) {
        mc->nb_vr_inner[j] = m[j];
        for (i = 0; i < m[j]; i++) {
            t = s[j][i];
            novc_eval_vokn(ts.inner[j], t, &temp);
            mc->vertex[k].u = temp.absi;
            mc->vertex[k].v = temp.ordo;
            mc->zt[k] = t;
            mc->flag[k] = j;
            k++;
        }
    }
    for (i = 0; i < nin; i++)
        free(s[i]);
    free(s);
    free(m);
    free(out);
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

