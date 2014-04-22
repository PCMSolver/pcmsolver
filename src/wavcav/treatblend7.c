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
#include "partcas.h"
#include "meshsas.h"


int trk = 0;



int pafn_read_wopc(side_flag SF, int z)
{
    if (z == ON_ALPHA)
        return SF.aflg;
    else if (z == ON_BETA)
        return SF.bflg;
    else if (z == ON_GAMMA)
        return SF.gflg;
    else if (z == ON_DELTA)
        return SF.dflg;
    else {
        fprintf(tmpout, "None of the four cases\n");
        exit(0);
    }
    return -1;
}



void werc_writ_vufh(side_flag * SF, int z, int value)
{
    if (z == ON_ALPHA)
        SF->aflg = value;
    else if (z == ON_BETA)
        SF->bflg = value;
    else if (z == ON_GAMMA)
        SF->gflg = value;
    else if (z == ON_DELTA)
        SF->dflg = value;
    else {
        fprintf(tmpout, "None of the four cases\n");
        exit(0);
    }
}



int hogm_find_givl_kocf(side_flag * SF, trmsrf * surf2, trmsrf * surf3, int w, int sd, prat_main_blend G, bd_box3D * B, int *exc, int *SD, double eps_inc)
{
    int nnd, ts_interf, z_w, z_i;
    int i, side, nd, found_id;
    double sml, dis;
    c_arc3D C;

    nnd = G.n_grs;
    z_w = G.N[w].trim_idx;
    if (G.N[w].type_node == T_TYPE) {
        fprintf(tmpout, "This applies only to Q_TYPE\n");
        exit(0);
    }
    if (sd == ON_ALPHA)
        poms_find_resk_lonb(surf2[z_w].pt.alpha, &C);
    if (sd == ON_BETA)
        poms_find_resk_lonb(surf2[z_w].pt.beta, &C);
    if (sd == ON_GAMMA)
        poms_find_resk_lonb(surf2[z_w].pt.gamma, &C);
    if (sd == ON_DELTA)
        poms_find_resk_lonb(surf2[z_w].pt.delta, &C);

    sml = LARGE_NUMBER;
    nd = -1;
    for (i = 0; i < nnd; i++)
        if ((i != w) && (exc[i] == 0)) {
            ts_interf = lafc_boun_gusd(B[i], B[w]);
            if (ts_interf == 1) {
                z_i = G.N[i].trim_idx;
                if (G.N[i].type_node == Q_TYPE) {
                    dis = jufq_dist_dunt(SF[i], surf2[z_i].pt, C, &side, &found_id);

                }
                if (G.N[i].type_node == T_TYPE)
                    dis = wujt_dist_kaqt(SF[i], surf3[z_i].pc, C, &side, &found_id);
                if ((found_id == 1) && (dis < sml)) {
                    sml = dis;
                    nd = i;
                    *SD = side;
                }
            }
        }
    if (nd != -1) {
        if (sml > eps_inc)
            nd = -1;
    }
    return nd;
}


int inc_node_blend_graph_interface(side_flag * SF, trmsrf * surf2, trmsrf * surf3, int w, int sd, prat_main_blend G, bd_box3D * B, int *exc, int *SD, double eps_inc)
{
    int nnd, ts_interf, z_w, z_i;
    int i, side, nd, found_id;
    double sml, dis;
    c_arc3D C;

    nnd = G.n_grs;
    z_w = G.N[w].trim_idx;
    if (G.N[w].type_node == Q_TYPE) {
        fprintf(tmpout, "This applies only to T_TYPE\n");
        exit(0);
    }
    if (sd == ON_ALPHA)
        poms_find_resk_lonb(surf3[z_w].pc.alpha, &C);
    if (sd == ON_BETA)
        poms_find_resk_lonb(surf3[z_w].pc.beta, &C);
    if (sd == ON_GAMMA)
        poms_find_resk_lonb(surf3[z_w].pc.gamma, &C);

    sml = LARGE_NUMBER;
    nd = -1;
    for (i = 0; i < nnd; i++)
        if ((i != w) && (exc[i] == 0)) {
            ts_interf = lafc_boun_gusd(B[i], B[w]);
            if (ts_interf == 1) {
                z_i = G.N[i].trim_idx;
                if (G.N[i].type_node == Q_TYPE)
                    dis = jufq_dist_dunt(SF[i], surf2[z_i].pt, C, &side, &found_id);
                if (G.N[i].type_node == T_TYPE)
                    dis = wujt_dist_kaqt(SF[i], surf3[z_i].pc, C, &side, &found_id);
                if ((found_id == 1) && (dis < sml)) {
                    sml = dis;
                    nd = i;
                    *SD = side;
                }
            }
        }
    if (nd != -1) {
        if (sml > eps_inc)
            nd = -1;
    }
    return nd;
}




int wumq_find_gopc_cunw(side_flag * SF, trmsrf * surf2, trmsrf * surf3, int w, int sd, prat_main_blend G, bd_box3D * B, int *exc, int *SD, double eps_inc)
{
    int nnd, ts_interf, z_w, z_i, i, side, nd, found_id;
    double sml, dis;
    c_arc3D C;

    nnd = G.n_grs;
    z_w = G.N[w].trim_idx;
    if (G.N[w].type_node == Q_TYPE) {
        fprintf(tmpout, "This applies only to T_TYPE\n");
        exit(0);
    }
    if (sd == ON_ALPHA)
        poms_find_resk_lonb(surf3[z_w].pc.alpha, &C);
    if (sd == ON_BETA)
        poms_find_resk_lonb(surf3[z_w].pc.beta, &C);
    if (sd == ON_GAMMA)
        poms_find_resk_lonb(surf3[z_w].pc.gamma, &C);

    sml = LARGE_NUMBER;
    nd = -1;
    for (i = 0; i < nnd; i++) {
        if ((i != w) && (exc[i] == 0) && (G.N[i].type_node == Q_TYPE)) {
            ts_interf = lafc_boun_gusd(B[i], B[w]);
            if (ts_interf == 1) {
                z_i = G.N[i].trim_idx;
                if (G.N[i].type_node == Q_TYPE)
                    dis = jufq_dist_dunt(SF[i], surf2[z_i].pt, C, &side, &found_id);
                if (G.N[i].type_node == T_TYPE)
                    dis = wujt_dist_kaqt(SF[i], surf3[z_i].pc, C, &side, &found_id);
                if (dis < sml) {
                    sml = dis;
                    nd = i;
                    *SD = side;
                }
            }
        }
    }
    if (nd != -1) {
        if (sml > eps_inc)
            nd = -1;
    }
    return nd;
}


int cigh_amon_zuvf(prat_main_blend G, int a, int b, int N)
{
    int i, n1, n2;
    for (i = 0; i < N; i++) {
        n1 = G.E[i].frvrt;
        n2 = G.E[i].scvrt;
        if ((n1 == a) && (n2 == b))
            return 1;
        if ((n1 == b) && (n2 == a))
            return 1;
    }
    return 0;
}


int cekz_test_hiqp(int *exc, prat_main_blend G)
{
    int nnd, w, res = 1;
    nnd = G.n_grs;
    for (w = 0; w < nnd; w++) {
        if ((G.N[w].val < 3) && (G.N[w].type_node == T_TYPE)) {
            res = 0;
            break;
        }
        if ((G.N[w].val < 2) && (G.N[w].type_node == Q_TYPE)) {
            res = 0;
            break;
        }
    }
    return res;
}



void jevm_trea_dufl(side_flag * SF, trmsrf * surf2, trmsrf * surf3, bd_box3D * B, prat_main_blend * G, int *exc, included_sides * inc, double eps_inc)
{
    int i, j, w, e, *sid, nd, val, fl;
    int ts_amo, q, n_sid, nnd, SD;
    e = G->k_grs;
    nnd = G->n_grs;
    sid = (int *) malloc(3 * sizeof(int));
    for (w = 0; w < nnd; w++)
        if ((exc[w] == 0) && (G->N[w].type_node == Q_TYPE)) {
            q = G->N[w].trim_idx;
            n_sid = 2;
            sid[0] = inc[w].side_1;
            sid[1] = inc[w].side_2;
            for (j = 0; j < n_sid; j++) {
                fl = pafn_read_wopc(SF[w], sid[j]);
                if (fl == -1) {
                    nd = hogm_find_givl_kocf(SF, surf2, surf3, w, sid[j], *G, B, exc, &SD, eps_inc);
                    if (nd != -1) {
                        ts_amo = cigh_amon_zuvf(*G, w, nd, e);
                        if (ts_amo == 0) {
                            G->E[e].frvrt = w;
                            G->E[e].scvrt = nd;
                            G->E[e].side_1 = sid[j];
                            G->E[e].side_2 = SD;
                            werc_writ_vufh(&SF[w], sid[j], 0);
                            werc_writ_vufh(&SF[nd], SD, 0);

                            val = G->N[w].val;
                            if (val >= 2) {
                                for (i = 0; i < 3; i++)
                                    fprintf(tmpout, "1-node[%d]  incident G-edge[%d]=%d\n", w, i, G->N[w].inc_edge[i]);
                                fprintf(tmpout, "Valence more than 3\n");
                                exit(0);
                            }
                            G->N[w].inc_edge[val] = e;
                            G->N[w].val = val + 1;

                            val = G->N[nd].val;
                            if (((G->N[nd].type_node == T_TYPE) && (val >= 3)) || ((G->N[nd].type_node == Q_TYPE) && (val >= 2))) {
                                for (i = 0; i < 3; i++)
                                    fprintf(tmpout, "2-node[%d]  indident G-edge[%d]=%d\n", nd, i, G->N[nd].inc_edge[i]);
                                fprintf(tmpout, "Valence more than 3\n");
                                exit(0);
                            }
                            G->N[nd].inc_edge[val] = e;
                            G->N[nd].val = val + 1;

                            if (G->N[w].val == 2)
                                exc[w] = +1;
                            if ((G->N[nd].val == 3) && (G->N[nd].type_node == T_TYPE))
                                exc[nd] = +1;
                            if ((G->N[nd].val == 2) && (G->N[nd].type_node == Q_TYPE))
                                exc[nd] = +1;
                            e++;
                        }
                    }
                }
            }
        }
    free(sid);
    G->k_grs = e;
}



void lunh_trea_jidh(side_flag * SF, trmsrf * surf2, trmsrf * surf3, bd_box3D * B, prat_main_blend * G, int *exc, included_sides * inc, double eps_inc)
{
    int i, j, w, e, *sid, nd, val, fl;
    int ts_amo, q, n_sid, nnd, SD;
    e = G->k_grs;
    nnd = G->n_grs;
    sid = (int *) malloc(3 * sizeof(int));
    for (w = 0; w < nnd; w++)
        if ((exc[w] == 0) && (G->N[w].type_node == T_TYPE)) {
            q = G->N[w].trim_idx;
            n_sid = 3;
            sid[0] = ON_ALPHA;
            sid[1] = ON_BETA;
            sid[2] = ON_GAMMA;
            for (j = 0; j < n_sid; j++) {
                fl = pafn_read_wopc(SF[w], sid[j]);
                if (fl == -1) {
                    nd = wumq_find_gopc_cunw(SF, surf2, surf3, w, sid[j], *G, B, exc, &SD, eps_inc);
                    if (nd != -1) {
                        ts_amo = cigh_amon_zuvf(*G, w, nd, e);
                        if (ts_amo == 0) {
                            G->E[e].frvrt = w;
                            G->E[e].scvrt = nd;
                            G->E[e].side_1 = sid[j];
                            G->E[e].side_2 = SD;
                            werc_writ_vufh(&SF[w], sid[j], 0);
                            werc_writ_vufh(&SF[nd], SD, 0);

                            val = G->N[w].val;
                            if (val >= 3) {
                                for (i = 0; i < 3; i++)
                                    fprintf(tmpout, "1-node[%d]  incident G-edge[%d]=%d\n", w, i, G->N[w].inc_edge[i]);
                                fprintf(tmpout, "Valence more than 3\n");
                                exit(0);
                            }
                            G->N[w].inc_edge[val] = e;
                            G->N[w].val = val + 1;

                            val = G->N[nd].val;
                            if (((G->N[nd].type_node == T_TYPE) && (val >= 3)) || ((G->N[nd].type_node == Q_TYPE) && (val >= 2))) {
                                for (i = 0; i < 3; i++)
                                    fprintf(tmpout, "2-node[%d]  indident G-edge[%d]=%d\n", nd, i, G->N[nd].inc_edge[i]);
                                fprintf(tmpout, "Valence more than 3\n");
                                exit(0);
                            }
                            G->N[nd].inc_edge[val] = e;
                            G->N[nd].val = val + 1;

                            if (G->N[w].val == 3)
                                exc[w] = +1;
                            if ((G->N[nd].val == 3) && (G->N[nd].type_node == T_TYPE))
                                exc[nd] = +1;
                            if ((G->N[nd].val == 2) && (G->N[nd].type_node == Q_TYPE))
                                exc[nd] = +1;
                            e++;
                        }
                    }
                }
            }
        }
    G->k_grs = e;
    free(sid);
}


void display_loc_sit(trmsrf surf, int sid)
{
    if (sid == ON_ALPHA)
        nepf_disp_bulp(surf.pc.alpha);
    if (sid == ON_BETA)
        nepf_disp_bulp(surf.pc.beta);
    if (sid == ON_GAMMA)
        nepf_disp_bulp(surf.pc.gamma);
}



void interface_inc_dead(int *wr, int *sd, int n_dd, int *map, side_flag * SF, trmsrf * surf2, trmsrf * surf3, bd_box3D * B, prat_main_blend * G, int *exc, included_sides * inc, double eps_inc)
{
    int i, w, e, sid, nd, val, fl;
    int ts_amo, q, nnd, SD, p, z;
    e = G->k_grs;
    nnd = G->n_grs;
    for (p = 0; p < n_dd; p++) {
        z = wr[p];
        w = map[z];
        q = G->N[w].trim_idx;
        sid = sd[p];


        fl = pafn_read_wopc(SF[w], sid);
        if (fl == -1) {
            trk = 1;
            nd = inc_node_blend_graph_interface(SF, surf2, surf3, w, sid, *G, B, exc, &SD, eps_inc);

            if (nd != -1) {
                ts_amo = cigh_amon_zuvf(*G, w, nd, e);
                if (ts_amo == 0) {
                    G->E[e].frvrt = w;
                    G->E[e].scvrt = nd;
                    G->E[e].side_1 = sid;
                    G->E[e].side_2 = SD;
                    fprintf(tmpout, "New interface edge[%d]=[%d,%d]\n", e, G->E[e].frvrt, G->E[e].scvrt);
                    werc_writ_vufh(&SF[w], sid, 0);
                    werc_writ_vufh(&SF[nd], SD, 0);

                    val = G->N[w].val;
                    if (val >= 3) {
                        for (i = 0; i < 3; i++)
                            fprintf(tmpout, "1-node[%d]  incident G-edge[%d]=%d\n", w, i, G->N[w].inc_edge[i]);
                        fprintf(tmpout, "Valence more than 3\n");
                        exit(0);
                    }
                    G->N[w].inc_edge[val] = e;
                    G->N[w].val = val + 1;

                    val = G->N[nd].val;
                    if (((G->N[nd].type_node == T_TYPE) && (val >= 3)) || ((G->N[nd].type_node == Q_TYPE) && (val >= 2))) {
                        for (i = 0; i < 3; i++)
                            fprintf(tmpout, "2-node[%d]  indident G-edge[%d]=%d\n", nd, i, G->N[nd].inc_edge[i]);
                        fprintf(tmpout, "Valence more than 3\n");
                        exit(0);
                    }
                    G->N[nd].inc_edge[val] = e;
                    G->N[nd].val = val + 1;

                    if (G->N[w].val == 3)
                        exc[w] = +1;
                    if ((G->N[nd].val == 3) && (G->N[nd].type_node == T_TYPE))
                        exc[nd] = +1;
                    if ((G->N[nd].val == 2) && (G->N[nd].type_node == Q_TYPE))
                        exc[nd] = +1;
                    e++;
                }
            }
        }
    }
    G->k_grs = e;
}


void display_side(int sd)
{
    if (sd == ON_ALPHA)
        fprintf(tmpout, "ON_ALPHA");
    if (sd == ON_BETA)
        fprintf(tmpout, "ON_BETA");
    if (sd == ON_GAMMA)
        fprintf(tmpout, "ON_GAMMA");
}



void loc_ent_include(trmsrf * surf3, side_flag * SF, hash_entry * gap, int n_gap, prat_main_blend * g_loc, prat_main_blend * G, int *map, int *exc)
{
    int p, ned_loc, i, j, n1, n2, k_ed, N1, N2, val;
    int w1, w2, nnd_loc, *map_e, w, nd, e, sd1, sd2;
    k_ed = 0;
    for (p = 0; p < n_gap; p++) {
        ned_loc = g_loc[p].k_grs;
        map_e = (int *) malloc(ned_loc * sizeof(int));

        for (i = 0; i < ned_loc; i++) {
            n1 = g_loc[p].E[i].frvrt;
            n2 = g_loc[p].E[i].scvrt;
            w1 = gap[p].list[n1];
            w2 = gap[p].list[n2];
            N1 = map[w1];
            N2 = map[w2];
            sd1 = g_loc[p].E[i].side_1;
            sd2 = g_loc[p].E[i].side_2;
            G->E[k_ed].frvrt = N1;
            G->E[k_ed].scvrt = N2;
            fprintf(tmpout, "New local edge[%d]=[%d,%d]\n", k_ed, N1, N2);
            if (0) {
                matc_disp_huqm(surf3[G->N[N1].trim_idx]);
                fprintf(tmpout, "========================\n");
                matc_disp_huqm(surf3[G->N[N2].trim_idx]);
                fprintf(tmpout, "========================\n");
                fprintf(tmpout, "sd1=");
                display_side(sd1);
                fprintf(tmpout, "\n");
                fprintf(tmpout, "sd2=");
                display_side(sd2);
                fprintf(tmpout, "\n");
            }

            G->E[k_ed].side_1 = sd1;
            G->E[k_ed].side_2 = sd2;
            werc_writ_vufh(&SF[N1], sd1, 0);
            werc_writ_vufh(&SF[N2], sd2, 0);
            map_e[i] = k_ed;
            k_ed++;
        }

        nnd_loc = g_loc[p].n_grs;
        for (i = 0; i < nnd_loc; i++) {
            w = gap[p].list[i];
            nd = map[w];
            val = g_loc[p].N[i].val;
            for (j = 0; j < val; j++) {
                e = g_loc[p].N[i].inc_edge[j];
                G->N[nd].inc_edge[j] = map_e[e];
            }
            G->N[nd].val = val;
            if (val == 3)
                exc[nd] = 1;
        }
        free(map_e);
    }
    G->k_grs = k_ed;
}


void det_edges_inc_simple(trmsrf * surf2, trmsrf * surf3, bd_box3D * B, prat_main_blend * G, int *exc, included_sides * inc, double eps_min, double eps_max, int n_trials)
{
    int i, ts, nnd;
    double eps_inc, step, lambda;
    side_flag *SF;
    nnd = G->n_grs;

    SF = (side_flag *) malloc(nnd * sizeof(side_flag));
    for (i = 0; i < nnd; i++) {
        SF[i].aflg = -1;
        SF[i].bflg = -1;
        SF[i].gflg = -1;
        SF[i].dflg = -1;
    }

    step = 1.0 / (double) n_trials;
    for (i = 0; i <= n_trials; i++) {
        fprintf(tmpout, "trial=%d\n", i);
        lambda = (double) i *step;
        eps_inc = lambda * eps_max + (1.0 - lambda) * eps_min;
        lunh_trea_jidh(SF, surf2, surf3, B, G, exc, inc, eps_inc);
        jevm_trea_dufl(SF, surf2, surf3, B, G, exc, inc, eps_inc);
        ts = cekz_test_hiqp(exc, *G);
        if (ts == 1)
            break;
    }
    free(SF);
}



void det_graph_bl_simple(sphere * S, trmsrf * surf2, int nb_surf2, trmsrf * surf3, int nb_surf3, blend_cpx BC, prat_main_blend * G, included_sides * inc)
{
    int i, j, k, N, w, M, nnd, side_1, side_2, p, q, z, *exc, n_trials = 10, suc;
    int ned, n_op, qw[4] = { ON_ALPHA, ON_BETA, ON_GAMMA, ON_DELTA };
    double marg = 0.1, eps = 1.0e-5, eps_min = 1.0e-5, eps_max = 1.0e-3;
    bd_box3D *B;

    N = BC.bt_grs;
    if (nb_surf2 != N) {
        fprintf(tmpout, "Requirement is not met\n");
        fprintf(tmpout, "[nb_surf2,nb_surf3]=[%d,%d]\n", nb_surf2, nb_surf3);
        fprintf(tmpout, "Nb of blend tors=%d\n", N);
        exit(0);
    }
    M = N + nb_surf3;
    B = (bd_box3D *) malloc(M * sizeof(bd_box3D));

    k = 0;
    for (i = 0; i < N; i++) {
        suc = jish_inci_gelh(S, surf2, BC, i, &side_1, &side_2, eps);
        if (suc == SUCCESS) {
            w = BC.BT[i].trim_idx;
            G->N[k].trim_idx = w;
            G->N[k].type_node = Q_TYPE;
            letd_boun_noth(surf2[w].pt, marg, &B[k]);
            inc[k].side_1 = -1;
            inc[k].side_2 = -1;
            for (j = 0; j < 4; j++) {
                z = qw[j];
                if ((z != side_1) && (z != side_2)) {
                    q = inc[k].side_1;
                    if (q == -1)
                        inc[k].side_1 = z;
                    else
                        inc[k].side_2 = z;
                }
            }
            k++;
        }
    }

    for (i = 0; i < nb_surf3; i++) {
        G->N[k].trim_idx = i;
        G->N[k].type_node = T_TYPE;
        nojp_boun_zovq(surf3[i].pc, marg, &B[k]);
        k++;
    }
    nnd = k;
    G->n_grs = nnd;
    fprintf(tmpout, "Number of nodes of blend prat_main=%d\n", nnd);

    exc = (int *) malloc(nnd * sizeof(int));
    for (i = 0; i < nnd; i++) {
        exc[i] = 0;
        G->N[i].val = 0;
    }
    G->k_grs = 0;
    det_edges_inc_simple(surf2, surf3, B, G, exc, inc, eps_min, eps_max, n_trials);

    free(exc);
    free(B);
    fprintf(tmpout, "Number of edges of blend prat_main=%d\n", G->k_grs);

    ned = G->k_grs;
    for (p = 0; p < ned; p++) {
        n_op = nevl_oppo_cikj(*G, p, G->E[p].opp);
        G->E[p].nb_opp = n_op;
    }
    fprintf(tmpout, "Search for opposite sides is complete\n");

}


void cetm_find_cojq_vurg(trmsrf * surf2, trmsrf * surf3, bd_box3D * B, prat_main_blend * G, int *exc, included_sides * inc, double eps_min, double eps_max, int n_trials)
{
    int i, ts, nnd;
    double eps_inc, step, lambda;
    side_flag *SF;
    nnd = G->n_grs;

    SF = (side_flag *) malloc(nnd * sizeof(side_flag));
    for (i = 0; i < nnd; i++) {
        SF[i].aflg = -1;
        SF[i].bflg = -1;
        SF[i].gflg = -1;
        SF[i].dflg = -1;
    }

    step = 1.0 / (double) n_trials;
    for (i = 0; i <= n_trials; i++) {
        fprintf(tmpout, "trial=%d\n", i);
        lambda = (double) i *step;
        eps_inc = lambda * eps_max + (1.0 - lambda) * eps_min;
        lunh_trea_jidh(SF, surf2, surf3, B, G, exc, inc, eps_inc);
        jevm_trea_dufl(SF, surf2, surf3, B, G, exc, inc, eps_inc);
        ts = cekz_test_hiqp(exc, *G);
        if (ts == 1)
            break;
    }
    free(SF);
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

