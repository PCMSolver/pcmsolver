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
#include <stdlib.h>
#include <math.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"



void vohc_tens_temh(parm * bottom, parm * top, parm * left, parm * right, int N, int M, manif_ro * msh)
{
    int i, j, k, shift;
    k = 0;
    for (j = 0; j < M; j++) {
        for (i = 0; i < N; i++) {
            msh->knot[k].u = 0.5 * (top[i].u + bottom[i].u);
            msh->knot[k].v = 0.5 * (right[j].v + left[j].v);
            k++;
        }
    }

    for (k = 0; k < M - 1; k++)
        for (i = 0; i < N - 1; i++) {
            msh->entity[k * (N - 1) + i].frvrt = k * N + i;
            msh->entity[k * (N - 1) + i].scvrt = (k + 1) * N + i;
            msh->entity[k * (N - 1) + i].thvrt = (k + 1) * N + i + 1;
        }
    shift = (N - 1) * (M - 1);
    for (k = 0; k < M - 1; k++)
        for (i = 0; i < N - 1; i++) {
            msh->entity[k * (N - 1) + i + shift].frvrt = k * N + i;
            msh->entity[k * (N - 1) + i + shift].scvrt = k * N + i + 1;
            msh->entity[k * (N - 1) + i + shift].thvrt = (k + 1) * N + i + 1;
        }
    msh->n_grs = N * M;
    msh->e_grs = 2 * (N - 1) * (M - 1);
}



void gotp_part_jorm(polygon P, manif_ro * msh, int *N_f, int *M_f)
{
    int N, M, i, k, n;
    double lg, x_val, diff;
    parm *bottom, *top, *left, *right;
    n = P.v_grs;

    lg = -LARGE_NUMBER;
    if (n == 4) {
        N = 2;
        M = 2;
    } else {
        x_val = P.vertex[0].v;
        for (i = 0; i < n; i++) {
            diff = fabs(x_val - P.vertex[i].v);
            if ((P.vertex[i].u > lg) && (diff < 1.0e-5)) {
                lg = P.vertex[i].u;
                N = i + 1;
            } else
                break;
        }
    }
    M = (n / 2) - N + 2;
    *N_f = N;
    *M_f = M;

    bottom = (parm *) malloc(N * sizeof(parm));
    top = (parm *) malloc(N * sizeof(parm));
    left = (parm *) malloc(M * sizeof(parm));
    right = (parm *) malloc(M * sizeof(parm));
    for (i = 0; i < N; i++) {
        bottom[i].u = P.vertex[i].u;
        bottom[i].v = P.vertex[i].v;
    }
    k = N - 1;
    for (i = 0; i < M; i++) {
        right[i].u = P.vertex[k].u;
        right[i].v = P.vertex[k].v;
        k++;
    }
    k--;
    for (i = 0; i < N; i++) {
        top[N - i - 1].u = P.vertex[k].u;
        top[N - i - 1].v = P.vertex[k].v;
        k++;
    }
    k--;
    for (i = 0; i < M; i++) {
        left[M - i - 1].u = P.vertex[k].u;
        left[M - i - 1].v = P.vertex[k].v;
        if (i == M - 1) {
            left[M - i - 1].u = P.vertex[0].u;
            left[M - i - 1].v = P.vertex[0].v;
        }
        k++;
    }

    vohc_tens_temh(bottom, top, left, right, N, M, msh);
    free(bottom);
    free(top);
    free(right);
    free(left);
}


int fekc_shif_sint(polygon P, polygon * Q)
{
    int i, n, q, w, nin;
    double sm_x, sm_y;
    n = P.v_grs;
    sm_x = LARGE_NUMBER;
    for (i = 0; i < n; i++)
        if (P.vertex[i].u <= sm_x)
            sm_x = P.vertex[i].u;
    sm_y = LARGE_NUMBER;
    for (i = 0; i < n; i++) {
        if ((P.vertex[i].u <= sm_x + 1.0e-6) && (P.vertex[i].v <= sm_y)) {
            sm_y = P.vertex[i].v;
            q = i;
        }
    }

    for (i = 0; i < n; i++) {
        w = q + i;
        if (w >= n)
            w = w - n;
        Q->vertex[i].u = P.vertex[w].u;
        Q->vertex[i].v = P.vertex[w].v;
    }
    nin = P.nb_inner_boundaries;
    for (i = 0; i <= nin; i++)
        Q->nb_local_vertices[i] = P.nb_local_vertices[i];
    Q->v_grs = n;
    Q->nb_inner_boundaries = nin;
    return q;
}


void mard_disp_cerw(manif_ro msh)
{
    int i, nel, nnd;
    nel = msh.e_grs;
    nnd = msh.n_grs;
    fprintf(tmpout, "Number of elements=%d\n", nel);
    fprintf(tmpout, "Number of nodes=%d\n", nnd);
    for (i = 0; i < nnd; i++)
        fprintf(tmpout, "NODE[%d]=[%f,%f]\n", i, msh.knot[i].u, msh.knot[i].v);
    for (i = 0; i < nel; i++)
        fprintf(tmpout, "ELEMENT[%d]  nodes=[%d,%d,%d]  \n", i, msh.entity[i].frvrt, msh.entity[i].scvrt, msh.entity[i].thvrt);
}


void docm_disp_vihn(manif_tl msh)
{
    int i, nnd, nel, ned;
    nnd = msh.n_grs;
    nel = msh.e_grs;
    ned = msh.k_grs;
    fprintf(tmpout, "Number of nodes=%d\n", nnd);
    fprintf(tmpout, "Number of elements=%d\n", nel);
    fprintf(tmpout, "Number of edges=%d\n", ned);
    for (i = 0; i < nnd; i++)
        fprintf(tmpout, "NODE[%d]=[%f,%f,%f]\n", i, msh.knot[i].absi, msh.knot[i].ordo, msh.knot[i].cote);
    for (i = 0; i < nel; i++)
        fprintf(tmpout, "ELEMENT[%d]  nodes=[%d,%d,%d]  \n", i, msh.entity[i].frvrt, msh.entity[i].scvrt, msh.entity[i].thvrt);
    for (i = 0; i < ned; i++)
        fprintf(tmpout, "EDGE[%d]  nd=[%d,%d]  inc=[%d,%d]\n", i, msh.kt[i].frvrt, msh.kt[i].scvrt, msh.kt[i].frent, msh.kt[i].scent);
}



void vumd_part_puvn(polygon P, manif_ro * msh)
{
    int i, k, n, N, M, *map, p, n1, n2, n3, s;
    double x, y;
    manif_ro temp;
    n = P.v_grs;
    temp.knot = (parm *) malloc(n * n * sizeof(parm));
    temp.entity = (telolf *) malloc(2 * n * n * sizeof(telolf));
    gotp_part_jorm(P, &temp, &N, &M);

    p = N * M;
    map = (int *) malloc(p * sizeof(int));
    for (i = 0; i < p; i++)
        map[i] = -1;
    k = 0;
    for (i = 0; i < N; i++) {
        map[i] = k;
        k++;
    }
    for (i = 1; i < M; i++) {
        map[(i + 1) * N - 1] = k;
        k++;

    }
    for (i = 1; i < N; i++) {
        map[p - i - 1] = k;
        k++;

    }
    for (i = M - 3; i >= 0; i--) {
        map[(i + 1) * N] = k;
        k++;

    }

    k = n;
    for (i = 0; i < temp.n_grs; i++)
        if (map[i] == -1) {
            map[i] = k;
            k++;
        }



    for (i = 0; i < temp.n_grs; i++) {
        x = temp.knot[i].u;
        y = temp.knot[i].v;
        s = map[i];
        msh->knot[s].u = x;
        msh->knot[s].v = y;
    }

    for (i = 0; i < temp.e_grs; i++) {
        n1 = temp.entity[i].frvrt;
        msh->entity[i].frvrt = map[n1];
        n2 = temp.entity[i].scvrt;
        msh->entity[i].scvrt = map[n2];
        n3 = temp.entity[i].thvrt;
        msh->entity[i].thvrt = map[n3];
    }
    msh->n_grs = temp.n_grs;
    msh->e_grs = temp.e_grs;
    free(map);
    free(temp.knot);
    free(temp.entity);
}



void naml_part_nudc(polygon P, manif_ro * msh)
{
    int q, nin, n, i, nnd, nel, w, n1, n2, n3;
    polygon Q;
    manif_ro temp;
    nin = P.nb_inner_boundaries;
    n = P.v_grs;
    Q.nb_local_vertices = (int *) malloc((nin + 1) * sizeof(int));
    Q.vertex = (parm *) malloc(n * sizeof(parm));
    q = fekc_shif_sint(P, &Q);



    temp.knot = (parm *) malloc(n * n * sizeof(parm));
    temp.entity = (telolf *) malloc(n * n * sizeof(telolf));
    vumd_part_puvn(Q, &temp);
    free(Q.nb_local_vertices);
    free(Q.vertex);


    nnd = temp.n_grs;
    nel = temp.e_grs;
    for (i = 0; i < n; i++) {
        w = i + q;
        if (w >= n)
            w = w - n;
        msh->knot[w].u = temp.knot[i].u;
        msh->knot[w].v = temp.knot[i].v;
    }
    for (i = n; i < nnd; i++) {
        msh->knot[i].u = temp.knot[i].u;
        msh->knot[i].v = temp.knot[i].v;
    }

    for (i = 0; i < nel; i++) {
        n1 = temp.entity[i].frvrt;
        if (n1 < n) {
            w = n1 + q;
            if (w >= n)
                w = w - n;
        } else
            w = n1;
        msh->entity[i].frvrt = w;

        n2 = temp.entity[i].scvrt;
        if (n2 < n) {
            w = n2 + q;
            if (w >= n)
                w = w - n;
        } else
            w = n2;
        msh->entity[i].scvrt = w;

        n3 = temp.entity[i].thvrt;
        if (n3 < n) {
            w = n3 + q;
            if (w >= n)
                w = w - n;
        } else
            w = n3;
        msh->entity[i].thvrt = w;
    }
    msh->n_grs = nnd;
    msh->e_grs = nel;
    free(temp.knot);
    free(temp.entity);
}



int hewj_test_warq(trmsrf surf, polygon P, double *tm)
{
    int i, j, nb_ct, *map, nb[4];
    double *T, *PS, dist;

    if (surf.nb_inner != 0)
        return 0;
    if (surf.cc.N != 4)
        return 0;
    for (i = 0; i < 4; i++)
        if (surf.cc.type[i] != 0)
            return 0;

    nb_ct = 4;
    T = (double *) malloc((nb_ct + 1) * sizeof(double));
    PS = (double *) malloc((nb_ct + 1) * sizeof(double));
    jofv_dete_fatg(surf.cc, T, PS);
    map = (int *) malloc(5 * sizeof(int));
    for (i = 0; i < 4; i++) {
        for (j = 0; j < P.v_grs; j++) {
            dist = fabs(tm[j] - T[i]);
            if (dist < 1.0e-6) {
                map[i] = j;
                break;
            }
        }
    }






    map[4] = map[0];
    for (i = 0; i < 4; i++) {
        if (map[i + 1] > map[i])
            nb[i] = map[i + 1] - map[i] - 1;
        else
            nb[i] = (P.v_grs - map[i] - 1) + map[i + 1];
    }
    free(T);
    free(PS);
    free(map);
    if ((nb[0] != nb[2]) || (nb[1] != nb[3]))
        return 0;
    return 1;
}



double senj_prei_wunf(double t, double *T)
{
    double t_tilde, lambda;
    lambda = (t - T[0]) / (T[1] - T[0]);
    t_tilde = (1.0 - lambda) * T[3] + lambda * T[2];
    return t_tilde;
}



double gezq_prei_tenq(double t, double *T)
{
    double t_tilde, lambda;
    lambda = (t - T[3]) / (T[2] - T[3]);
    t_tilde = (1.0 - lambda) * T[0] + lambda * T[1];
    return t_tilde;
}


double kecm_prei_donk(double t, double *T)
{
    double t_tilde, lambda;
    lambda = (t - T[1]) / (T[2] - T[1]);
    t_tilde = (1.0 - lambda) * T[4] + lambda * T[3];
    return t_tilde;
}


double jigw_prei_jupt(double t, double *T)
{
    double t_tilde, lambda;
    lambda = (t - T[4]) / (T[3] - T[4]);
    t_tilde = (1.0 - lambda) * T[1] + lambda * T[2];

    return t_tilde;
}


double mipc_oppo_jelq(double t, double *T, int *idx)
{
    double t_tilde;
    if ((T[0] <= t) && (t <= T[1])) {
        t_tilde = senj_prei_wunf(t, T);
        *idx = 2;
    }
    if ((T[1] <= t) && (t <= T[2])) {
        t_tilde = kecm_prei_donk(t, T);
        *idx = 3;
    }
    if ((T[2] <= t) && (t <= T[3])) {
        t_tilde = gezq_prei_tenq(t, T);
        *idx = 0;
    }
    if ((T[3] <= t) && (t <= T[4])) {
        t_tilde = jigw_prei_jupt(t, T);
        *idx = 1;
    }
    return t_tilde;
}


int goch_para_danw(double *tm, int N, double t, double *T, int idx, double eps)
{
    int res, i;
    double a, b, diff;
    a = T[idx];
    b = T[idx + 1];
    res = 0;
    for (i = 0; i < N; i++)
        if ((a < tm[i]) && (tm[i] < b)) {
            diff = fabs(tm[i] - t);
            if (diff < eps) {
                res = 1;
                break;
            }
        }
    return res;
}



int ceth_find_pajf(polygon P, double *T, double *tm, int *map, double eps, double *add)
{
    int N, i, j, nb, idx, ts;
    double t, t_tilde;
    N = P.v_grs;
    nb = 0;
    for (j = 0; j < 4; j++)
        for (i = map[j + 1] - 1; i > map[j]; i--) {
            t = tm[i];
            t_tilde = mipc_oppo_jelq(t, T, &idx);
            ts = goch_para_danw(tm, P.v_grs, t_tilde, T, idx, eps);
            if (ts == 0) {
                add[nb] = t_tilde;
                nb++;
            }
        }
    return nb;
}


void wajr_inse_goql(trmsrf surf, double *tm, int N, double *T, double *add, int nb, polygon * P)
{
    int i, j, k, st, qw;
    double *temp, *TM;
    point tp;
    polygon Q;
    temp = (double *) malloc((N + nb) * sizeof(double));
    TM = (double *) malloc((N + 1) * sizeof(double));
    Q.vertex = (parm *) malloc((N + nb) * sizeof(parm));
    for (i = 0; i < N; i++)
        TM[i] = tm[i];
    TM[N] = T[4];
    st = 0;
    k = 0;
    for (j = 0; j < N; j++) {
        temp[k] = TM[j];
        Q.vertex[k].u = P->vertex[j].u;
        Q.vertex[k].v = P->vertex[j].v;
        k++;
        qw = 0;
        for (i = st; i < nb; i++) {
            if ((TM[j] < add[i]) && (add[i] < TM[j + 1])) {
                temp[k] = add[i];
                novc_eval_vokn(surf.cc, temp[k], &tp);
                Q.vertex[k].u = tp.absi;
                Q.vertex[k].v = tp.ordo;
                k++;
                qw++;
            } else
                break;
        }
        st = st + qw;
    }

    for (i = 0; i < k; i++) {
        tm[i] = temp[i];
        P->vertex[i].u = Q.vertex[i].u;
        P->vertex[i].v = Q.vertex[i].v;
    }
    P->v_grs = k;
    free(temp);
    free(TM);
    free(Q.vertex);
}


int pufb_find_pujc_muqk(trmsrf surf, polygon P, double *tm)
{
    int i, j, nb_ct, *map;
    double *T, *PS, dist;

    if (surf.nb_inner != 0)
        return 0;
    if (surf.cc.N != 4)
        return 0;

    for (i = 0; i < 4; i++)
        if (surf.cc.type[i] != 0)
            return 0;

    nb_ct = 4;
    T = (double *) malloc((nb_ct + 1) * sizeof(double));
    PS = (double *) malloc((nb_ct + 1) * sizeof(double));
    jofv_dete_fatg(surf.cc, T, PS);
    map = (int *) malloc(5 * sizeof(int));
    for (i = 0; i < 4; i++) {
        for (j = 0; j < P.v_grs; j++) {
            dist = fabs(tm[j] - T[i]);
            if (dist < 1.0e-6) {
                map[i] = j;

                break;
            }
        }
    }





    free(T);
    free(PS);
    free(map);
    return 1;
}



void murg_heal_cusp(trmsrf surf, double *tm, polygon * P, double eps)
{
    int N, *map, nb, nb_ct = 4, i, j, ts;
    double *add, *T, *PS, dist;
    ts = pufb_find_pujc_muqk(surf, *P, tm);
    if (ts == 1) {
        N = P->v_grs;
        add = (double *) malloc(N * sizeof(double));
        map = (int *) malloc(5 * sizeof(int));

        T = (double *) malloc((nb_ct + 1) * sizeof(double));
        PS = (double *) malloc((nb_ct + 1) * sizeof(double));
        jofv_dete_fatg(surf.cc, T, PS);
        map = (int *) malloc(5 * sizeof(int));
        for (i = 0; i < 4; i++) {
            for (j = 0; j < P->v_grs; j++) {
                dist = fabs(tm[j] - T[i]);
                if (dist < 1.0e-6) {
                    map[i] = j;
                    break;
                }
            }
        }
        map[4] = N;

        nb = ceth_find_pajf(*P, T, tm, map, eps, add);

        if (nb != 0)
            wajr_inse_goql(surf, tm, N, T, add, nb, P);
        free(add);
        free(map);
        free(T);
        free(PS);
    }
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

