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
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"
#include "sas.h"
#include "splinemol.h"
#include "meshsas.h"



void disc_parm_supporting(sphere * S, trmsrf ts, trmsrf * surf2, blend_cpx BC, int z, set_arcs SA, int *n_len, int *dis_ex, int **dis_in, int *endp, arc_supp A, supporting * sp)
{
    int N, i, j, k, w, *p_sa, val, q;
    int nx, n_comp, nin;
    int *exc, k_nw;
    double dis, sml;
    point *ss3D, X, Y;

    N = SA.ar_grs;
    val = BC.HE[z].nb;
    p_sa = (int *) malloc(N * sizeof(int));
    for (i = 0; i < N; i++) {
        w = A.q_sa[i];
        p_sa[i] = n_len[w];
    }

    exc = (int *) malloc(N * sizeof(int));
    for (i = 0; i < N; i++)
        exc[i] = 0;
    n_comp = ts.cc.N;
    if (n_comp == 2)
        fprintf(tmpout, "WARNING: cc has only two components\n");
    ss3D = (point *) malloc(n_comp * sizeof(point));
    sofl_segm_salc(ts, ts.cc, ss3D);
    k_nw = 0;
    for (i = 0; i < n_comp; i++) {
        nx = i + 1;
        if (nx == n_comp)
            nx = 0;
        getf_find_rogc_todj(ss3D[i], &X);
        getf_find_rogc_todj(ss3D[nx], &Y);
        sml = LARGE_NUMBER;
        q = -1;
        for (j = 0; j < N; j++)
            if (exc[j] == 0) {
                dis = cutj_dist_rulb(X, Y, SA.C[j]);
                if (dis < sml) {
                    sml = dis;
                    q = j;
                }
            }
        if (q == -1) {
            fprintf(tmpout, "Warning: unable to find bounding[X,Y]\n");
            exit(0);
        }
        sp->supp_ext[i] = q;

        dis_ex[i] = p_sa[q];
        exc[q] = 1;
        endp[k_nw] = A.q_sa[q];
        k_nw++;
    }

    nin = ts.nb_inner;
    for (k = 0; k < nin; k++) {
        n_comp = ts.inner[k].N;
        ss3D = (point *) malloc(n_comp * sizeof(point));
        sofl_segm_salc(ts, ts.inner[k], ss3D);
        for (i = 0; i < n_comp; i++) {
            nx = i + 1;
            if (nx == n_comp)
                nx = 0;
            getf_find_rogc_todj(ss3D[i], &X);
            getf_find_rogc_todj(ss3D[nx], &Y);
            sml = LARGE_NUMBER;
            q = -1;
            for (j = 0; j < N; j++)
                if (exc[j] == 0) {
                    dis = cutj_dist_rulb(X, Y, SA.C[j]);
                    if (dis < sml) {
                        sml = dis;
                        q = j;
                    }
                }
            if (q == -1) {
                fprintf(tmpout, "Warning: unable to find arc\n");
                exit(0);
            }
            sp->supp_int[k][i] = q;
            dis_in[k][i] = p_sa[q];
            exc[q] = 1;
            endp[k_nw] = A.q_sa[q];
            k_nw++;
        }
    }
    free(exc);
    free(p_sa);
    free(ss3D);
    sp->nin = ts.nb_inner;
}




void gen_loc_mult_conn(atom * S, int nb_sph, blend_cpx BC, set_arcs * SA, trmsrf * surf1, int q, trmsrf * surf2, int *supp, int *n_len, arc_supp * A, mult_conn * P, add_mc * P_app, supporting * sp, int max_nvt)
{
    int j, z, nc, nb_cr, nb_loc, ts_supp;
    int nin, *dis_ex, **dis_in, *endp;
    double acc;

    ts_supp = salr_test_jofl(surf1[q], S[supp[q]], 1.0e-4, &acc);
    if (ts_supp == 0) {
        fprintf(tmpout, "Not a supporting atom\n");
        exit(0);
    }

    nc = surf1[q].cc.N;
    nb_loc = nc;
    dis_ex = (int *) malloc(nc * sizeof(int));
    nin = surf1[q].nb_inner;
    dis_in = (int **) malloc(nin * sizeof(int *));
    for (j = 0; j < nin; j++) {
        nc = surf1[q].inner[j].N;
        dis_in[j] = (int *) malloc(nc * sizeof(int));
        nb_loc = nb_loc + nc;
    }
    z = supp[q];
    endp = (int *) malloc(nb_loc * sizeof(int));

    disc_parm_supporting(S, surf1[q], surf2, BC, z, SA[z], n_len, dis_ex, dis_in, endp, A[z], sp);
    free(endp);

    nb_cr = mult_conn_simple(surf1[q], dis_ex, dis_in, P, P_app, max_nvt);
    for (j = 0; j < nin; j++)
        free(dis_in[j]);
    free(dis_in);
    free(dis_ex);
}



int edges_interference(mult_conn mc, add_mc AM, int *list, int *forc_term)
{
    int n_ed1, n_ed2, n_ed, ts, dummy, f_trm;
    int *list1, *list2, ned, i;

    *forc_term = 0;
    ned = AM.nb_edge_mc;
    list1 = (int *) malloc(ned * sizeof(int));
    list2 = (int *) malloc(ned * sizeof(int));
    n_ed1 = edges_self_inters(mc, AM, list1);
    if (n_ed1 >= 1)
        fprintf(tmpout, "self-intersection detected\n");
    n_ed2 = 0;
    if (n_ed1 == 0) {
        n_ed2 = edges_coarse_ext(mc, AM, list2, &f_trm);
        if (f_trm == 1) {
            fprintf(tmpout, "force term: edges_coarse_ext() in edges_interference()\n");
            *forc_term = 1;
            free(list1);
            free(list2);
            return 0;
        }
    }

    n_ed = 0;
    for (i = 0; i < n_ed1; i++) {
        ts = gonl_arra_govj(list, n_ed, list1[i], &dummy);
        if (ts == 0) {
            list[n_ed] = list1[i];
            n_ed++;
        }
    }

    for (i = 0; i < n_ed2; i++) {
        ts = gonl_arra_govj(list, n_ed, list2[i], &dummy);
        if (ts == 0) {
            list[n_ed] = list2[i];
            n_ed++;
        }
    }

    free(list1);
    free(list2);
    return n_ed;
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

