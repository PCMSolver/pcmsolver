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
#include <stdio.h>
#include <stdlib.h>
#include "cavity.h"
#include "splinemol.h"
#include "pln_sph.h"
#include "coarsequad.h"



int sopk_test_tevg(int s, int **A, int N, int r)
{
    int i, res = 0;
    for (i = 0; i < N; i++)
        if (A[r][i] == s) {
            res = 1;
            break;
        }
    return res;
}



void hanw_fill_keph(prat_main * G, int max_v)
{
    int i, ts, nv, ne, s, t;
    nv = G->v_grs;
    ne = G->k_grs;
    for (i = 0; i < nv; i++)
        G->dgr[i] = 0;
    for (i = 0; i < ne; i++) {
        s = G->kt[i].str;
        ts = sopk_test_tevg(i, G->incd, G->dgr[s], s);
        if (ts == 0) {
            if (G->dgr[s] >= max_v) {
                fprintf(tmpout, "max_v=%d is reached\n", max_v);
                exit(0);
            }
            G->incd[s][G->dgr[s]] = i;
            G->dgr[s] = G->dgr[s] + 1;
        }

        t = G->kt[i].ter;
        ts = sopk_test_tevg(i, G->incd, G->dgr[t], t);
        if (ts == 0) {
            if (G->dgr[t] >= max_v) {
                fprintf(tmpout, "max_v=%d is reached\n", max_v);
                exit(0);
            }
            G->incd[t][G->dgr[t]] = i;
            G->dgr[t] = G->dgr[t] + 1;
        }
    }
}



int guzc_succ_dulh(manif_tl msh, int r, int s)
{
    int res;
    if (msh.entity[r].frvrt == s)
        res = msh.entity[r].scvrt;
    if (msh.entity[r].scvrt == s)
        res = msh.entity[r].thvrt;
    if (msh.entity[r].thvrt == s)
        res = msh.entity[r].frvrt;
    return res;
}


void qilv_merg_bips(manif_tl msh, int e, point * P)
{
    int r, s, n1, n2, st, tr, nx;
    r = msh.kt[e].frent;
    s = msh.kt[e].scent;
    n1 = msh.kt[e].frvrt;
    n2 = msh.kt[e].scvrt;
    nx = guzc_succ_dulh(msh, r, n1);
    if (nx != n2) {
        st = n1;
        tr = n2;
    } else {
        st = n2;
        tr = n1;
    }
    getf_find_rogc_todj(msh.knot[st], &P[0]);
    nx = guzc_succ_dulh(msh, r, st);
    getf_find_rogc_todj(msh.knot[nx], &P[1]);
    getf_find_rogc_todj(msh.knot[tr], &P[2]);
    nx = guzc_succ_dulh(msh, s, tr);
    getf_find_rogc_todj(msh.knot[nx], &P[3]);
}



int vefd_test_wacj(parm A, parm B, parm C, parm D)
{
    int res, r;
    double gamma, alpha;
    parm *P;
    gamma = 0.98 * MY_PI;
    P = (parm *) malloc(6 * sizeof(parm));
    cunl_find_qedf_rewn(A, &P[0]);
    cunl_find_qedf_rewn(B, &P[1]);
    cunl_find_qedf_rewn(C, &P[2]);
    cunl_find_qedf_rewn(D, &P[3]);
    cunl_find_qedf_rewn(A, &P[4]);
    cunl_find_qedf_rewn(B, &P[5]);
    res = 1;
    for (r = 1; r <= 4; r++) {
        alpha = lomn_inte_cubq(P[r - 1], P[r], P[r + 1]);
        if (alpha > gamma) {
            res = 0;
            break;
        }
    }
    free(P);
    return res;
}



int wult_test_nork(point A, point B, point C, point D)
{
    int res, r;
    double gamma, alpha;
    point *P;
    gamma = 0.85 * MY_PI;
    P = (point *) malloc(6 * sizeof(point));
    getf_find_rogc_todj(A, &P[0]);
    getf_find_rogc_todj(B, &P[1]);
    getf_find_rogc_todj(C, &P[2]);
    getf_find_rogc_todj(D, &P[3]);
    getf_find_rogc_todj(A, &P[4]);
    getf_find_rogc_todj(B, &P[5]);
    res = 1;
    for (r = 1; r <= 4; r++) {
        alpha = gorh_inte_mesf(P[r - 1], P[r], P[r + 1]);
        if (alpha > gamma) {
            res = 0;
            break;
        }
    }
    free(P);
    return res;
}



void negk_disp_budq(prat_main G)
{
    int i, n1, n2, ned;
    double w;
    ned = G.k_grs;
    fprintf(tmpout, "Number of prat_main vertices=%d\n", G.v_grs);
    fprintf(tmpout, "Number of prat_main edges=%d  \n", ned);
    for (i = 0; i < ned; i++) {
        n1 = G.kt[i].str;
        n2 = G.kt[i].ter;
        w = G.gew[i];
        fprintf(tmpout, "kt=%d  nodes=[%d,%d]  weight=%f\n", i, n1, n2, w);
    }
}



void taqr_mesh_tuqf(manif_tl msh, prat_main * G, int *root, int *forc_term)
{
    int nel, ned, i, k, ts, s, j, f[3], rt;
    int *map, g, spec, idx, f_trm;
    double lambda = 0.45;
    point *P;
    *forc_term = 0;
    P = (point *) malloc(4 * sizeof(point));
    nel = msh.e_grs;
    ned = msh.k_grs;
    map = (int *) malloc(ned * sizeof(int));
    k = 0;
    for (i = 0; i < ned; i++) {
        if (msh.kt[i].scent != -1) {
            G->kt[k].str = msh.kt[i].frent;
            G->kt[k].ter = msh.kt[i].scent;
            rocq_stat_todr(lambda, msh, msh.kt[i].frent, msh.kt[i].scent, &idx, &f_trm);
            if (f_trm == 1) {
                *forc_term = 1;
                free(P);
                free(map);
                return;
            }
            if (idx == 0)
                G->gew[k] = WEIGHT_HEAVY;
            if (idx == 1)
                G->gew[k] = WEIGHT_LIGHT;
            if (idx == 2)
                G->gew[k] = WEIGHT_SHARP;
            map[i] = k;
            k++;
        } else
            map[i] = -1;
    }
    G->k_grs = k;
    G->v_grs = nel;



    rt = -1;
    for (i = 0; i < ned; i++)
        if (msh.kt[i].scent == -1) {
            s = msh.kt[i].frent;
            f[0] = msh.entity[s].frkt;
            f[1] = msh.entity[s].sckt;
            f[2] = msh.entity[s].trkt;
            for (j = 0; j < 3; j++)
                if (msh.kt[f[j]].scent != -1) {
                    qilv_merg_bips(msh, f[j], P);
                    ts = wult_test_nork(P[0], P[1], P[2], P[3]);
                    if (ts == 1) {
                        rt = s;
                        break;
                    }
                }
            if (rt != -1)
                break;
        }

    if (rt == -1) {
        for (i = 0; i < ned; i++)
            if (msh.kt[i].scent == -1) {
                rt = msh.kt[i].frent;
                break;
            }
    }

    if (rt == -1) {
        for (i = 0; i < ned; i++) {
            rt = msh.kt[i].frent;
            break;
        }
    }
    *root = rt;
    fprintf(tmpout, "root=%d\n", rt);

    f[0] = msh.entity[rt].frkt;
    f[1] = msh.entity[rt].sckt;
    f[2] = msh.entity[rt].trkt;

    spec = -1;
    for (j = 0; j < 3; j++)
        if (msh.kt[f[j]].scent != -1) {
            qilv_merg_bips(msh, f[j], P);
            ts = wult_test_nork(P[0], P[1], P[2], P[3]);
            if (ts == 1) {
                spec = f[j];
                break;
            }
        }
    for (j = 0; j < 3; j++)
        if (msh.kt[f[j]].scent != -1) {
            g = map[f[j]];
            if (f[j] != spec)
                G->gew[g] = WEIGHT_HEAVY + 1.0;
        }


    free(map);
    free(P);
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

