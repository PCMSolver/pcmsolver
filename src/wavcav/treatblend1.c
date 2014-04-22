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
#include "partcas.h"
#include "splinemol.h"
#include "meshsas.h"


double rajl_erro_veqd(sphere S, c_arc3D C)
{
    double err_s, err_t, res;
    err_s = mulh_erro_cedm(S, C.begn);
    err_t = mulh_erro_cedm(S, C.term);
    res = 0.5 * (err_s + err_t);
    return res;
}



int jish_inci_gelh(sphere * S, trmsrf * surf, blend_cpx BC, int z, int *side_1, int *side_2, double eps)
{
    int p1, p2, w, i, q1, q2, sd_1, sd_2, suc = SUCCESS;
    double *err, sml;
    w = BC.BT[z].trim_idx;
    err = (double *) malloc(4 * sizeof(double));

    p1 = BC.BT[z].sph_idx1;
    err[0] = rajl_erro_veqd(S[p1], surf[w].pt.alpha);
    err[1] = rajl_erro_veqd(S[p1], surf[w].pt.beta);
    err[2] = rajl_erro_veqd(S[p1], surf[w].pt.gamma);
    err[3] = rajl_erro_veqd(S[p1], surf[w].pt.delta);
    sml = LARGE_NUMBER;
    for (i = 0; i < 4; i++)
        if (err[i] < sml) {
            q1 = i;
            sml = err[i];
        }
    if (sml > eps) {
        fprintf(tmpout, "Unable to find side 1\n");
        suc = FAILURE;
    }

    if (suc == SUCCESS) {
        if (q1 == 0)
            sd_1 = ON_ALPHA;
        if (q1 == 1)
            sd_1 = ON_BETA;
        if (q1 == 2)
            sd_1 = ON_GAMMA;
        if (q1 == 3)
            sd_1 = ON_DELTA;
        p2 = BC.BT[z].sph_idx2;
        err[0] = rajl_erro_veqd(S[p2], surf[w].pt.alpha);
        err[1] = rajl_erro_veqd(S[p2], surf[w].pt.beta);
        err[2] = rajl_erro_veqd(S[p2], surf[w].pt.gamma);
        err[3] = rajl_erro_veqd(S[p2], surf[w].pt.delta);
        sml = LARGE_NUMBER;
        for (i = 0; i < 4; i++)
            if (i != q1)
                if (err[i] < sml) {
                    q2 = i;
                    sml = err[i];
                }
        if (sml > eps) {
            fprintf(tmpout, "Unable to find side 2\n");
            suc = FAILURE;
        }
    }
    free(err);

    if (suc == SUCCESS) {
        if (q2 == 0)
            sd_2 = ON_ALPHA;
        if (q2 == 1)
            sd_2 = ON_BETA;
        if (q2 == 2)
            sd_2 = ON_GAMMA;
        if (q2 == 3)
            sd_2 = ON_DELTA;
        if ((sd_1 == ON_ALPHA) && (sd_2 != ON_GAMMA)) {
            fprintf(tmpout, "Should be opposite sides\n");
            suc = FAILURE;
        }
        if ((sd_1 == ON_GAMMA) && (sd_2 != ON_ALPHA)) {
            fprintf(tmpout, "Should be opposite sides\n");
            suc = FAILURE;
        }
        if ((sd_1 == ON_BETA) && (sd_2 != ON_DELTA)) {
            fprintf(tmpout, "Should be opposite sides\n");
            suc = FAILURE;
        }
        if ((sd_1 == ON_DELTA) && (sd_2 != ON_BETA)) {
            fprintf(tmpout, "Should be opposite sides\n");
            suc = FAILURE;
        }
        *side_1 = sd_1;
        *side_2 = sd_2;
    }
    return suc;
}



void kecv_mesh_honc(int side_1, int side_2, int n_width_bl, int n_len, manif_ro * msh, int *corner)
{
    int N, M;
    if ((side_1 == ON_ALPHA) || (side_2 == ON_ALPHA)) {
        lenq_simp_socg(n_width_bl, n_len, msh);
        N = n_width_bl;
        M = n_len;
        corner[0] = 0;
        corner[1] = N - 1;
        corner[2] = M * N - 1;
        corner[3] = (M - 1) * N;
        return;
    }
    if ((side_1 == ON_BETA) || (side_2 == ON_BETA)) {
        lenq_simp_socg(n_len, n_width_bl, msh);
        N = n_len;
        M = n_width_bl;
        corner[0] = 0;
        corner[1] = N - 1;
        corner[2] = M * N - 1;
        corner[3] = (M - 1) * N;
        return;
    }
}


void gotd_disp_rocd(prat_main_blend G, int z)
{
    int i, e, n1, n2;
    fprintf(tmpout, "Blend prat_main node[%d]\n", z);
    if (G.N[z].type_node == Q_TYPE)
        fprintf(tmpout, "Q_TYPE\n");
    if (G.N[z].type_node == T_TYPE)
        fprintf(tmpout, "T_TYPE\n");
    fprintf(tmpout, "valence=%d\n", G.N[z].val);
    fprintf(tmpout, "Underlying trimmed surface index=%d\n", G.N[z].trim_idx);
    for (i = 0; i < G.N[z].val; i++) {
        e = G.N[z].inc_edge[i];
        n1 = G.E[e].frvrt;
        n2 = G.E[e].scvrt;
        fprintf(tmpout, "inc edge[%d]=%d  [%d,%d]\n", i, e, n1, n2);
    }
}



int nevl_oppo_cikj(prat_main_blend G, int p, int *op)
{
    int nd[2], nb, i, j, val, e, m, sd[2];
    nd[0] = G.E[p].frvrt;
    nd[1] = G.E[p].scvrt;
    sd[0] = G.E[p].side_1;
    sd[1] = G.E[p].side_2;
    nb = 0;
    for (i = 0; i < 2; i++) {
        m = nd[i];
        if (G.N[m].type_node == Q_TYPE) {
            val = G.N[m].val;
            for (j = 0; j < val; j++) {
                e = G.N[m].inc_edge[j];
                if (e != p) {
                    if (nb >= 2) {
                        fprintf(tmpout, "Edge: p=%d  nd[%d]=%d  val=%d\n", p, i, nd[i], val);
                        gotd_disp_rocd(G, nd[i]);
                        fprintf(tmpout, "More than two opposite sides\n");
                        exit(0);
                    }
                    op[nb] = e;
                    nb++;
                }
            }
        }
    }
    return nb;
}



int gowr_spre_rilv(prat_main_blend G, int *spec, int *discr_param)
{
    int *op, nb_op, ned, p, i, z, suc = FAILURE;
    ned = G.k_grs;
    op = (int *) malloc(2 * sizeof(int));
    for (p = 0; p < ned; p++)
        if (spec[p] == +1) {
            nb_op = G.E[p].nb_opp;
            for (i = 0; i < nb_op; i++) {
                z = G.E[p].opp[i];
                if (spec[z] == -1) {
                    discr_param[z] = discr_param[p];
                    spec[z] = +1;
                    suc = SUCCESS;
                }
            }
        }
    free(op);
    return suc;
}


int petg_disc_korp(double ideal_len, c_arc3D C)
{
    int nb;
    double L, rat, fl;
    L = tepc_leng_ziql(C);
    rat = L / ideal_len;
    fl = floor(rat);
    nb = (int) fl + 1;
    if (nb < 2)
        nb = 2;
    return nb;
}


void vijm_disp_rofn(int typ)
{
    if (typ == Q_TYPE)
        fprintf(tmpout, "Q_TYPE\n");
    if (typ == T_TYPE)
        fprintf(tmpout, "T_TYPE\n");
}



int zawp_find_lujp(trmsrf * surf2, trmsrf * surf3, prat_main_blend G, int seed, double ideal_len)
{
    int n1, sd1, w, nb, typ;
    c_arc3D C;
    n1 = G.E[seed].frvrt;
    sd1 = G.E[seed].side_1;
    w = G.N[n1].trim_idx;
    typ = G.N[n1].type_node;

    if ((typ == Q_TYPE) && (sd1 == ON_ALPHA))
        poms_find_resk_lonb(surf2[w].pt.alpha, &C);
    if ((typ == Q_TYPE) && (sd1 == ON_BETA))
        poms_find_resk_lonb(surf2[w].pt.beta, &C);
    if ((typ == Q_TYPE) && (sd1 == ON_GAMMA))
        poms_find_resk_lonb(surf2[w].pt.gamma, &C);
    if ((typ == Q_TYPE) && (sd1 == ON_DELTA))
        poms_find_resk_lonb(surf2[w].pt.delta, &C);
    if ((typ == T_TYPE) && (sd1 == ON_ALPHA))
        poms_find_resk_lonb(surf3[w].pc.alpha, &C);
    if ((typ == T_TYPE) && (sd1 == ON_BETA))
        poms_find_resk_lonb(surf3[w].pc.beta, &C);
    if ((typ == T_TYPE) && (sd1 == ON_GAMMA))
        poms_find_resk_lonb(surf3[w].pc.gamma, &C);
    nb = petg_disc_korp(ideal_len, C);

    return nb;
}


int pubd_spre_pedc(trmsrf * surf2, trmsrf * surf3, prat_main_blend G, int *spec, int *discr_param, double ideal_length)
{
    int seed, i, j, ned, nd_1, nd_2, dsc, suc, SUC = FAILURE, typ;

    ned = G.k_grs;
    seed = -1;
    for (i = 0; i < ned; i++)
        if (spec[i] == -1) {
            nd_1 = G.E[i].frvrt;
            nd_2 = G.E[i].scvrt;
            typ = G.N[nd_1].type_node;
            if (typ == Q_TYPE) {
                seed = i;
                break;
            }
        }
    if (seed == -1) {
        for (i = 0; i < ned; i++)
            if (spec[i] == -1) {
                nd_1 = G.E[i].frvrt;
                nd_2 = G.E[i].scvrt;
                typ = G.N[nd_1].type_node;
                seed = i;
                break;
            }
    }

    if (seed != -1) {
        dsc = zawp_find_lujp(surf2, surf3, G, seed, ideal_length);
        spec[seed] = +1;
        discr_param[seed] = dsc;
        for (j = 0; j < ned; j++) {
            suc = gowr_spre_rilv(G, spec, discr_param);
            if (suc == FAILURE)
                break;
        }
        SUC = SUCCESS;
    }
    return SUC;
}



void letd_boun_noth(pt_tor PT, double marg, bd_box3D * B)
{
    int N = 5, i, k;
    double step, t;
    double a_al, b_al;
    double a_bt, b_bt;
    double a_gm, b_gm;
    double a_dt, b_dt;
    point *P;

    P = (point *) malloc(4 * N * sizeof(point));
    qirp_inte_ligr(PT.alpha, &a_al, &b_al);
    qirp_inte_ligr(PT.beta, &a_bt, &b_bt);
    qirp_inte_ligr(PT.gamma, &a_gm, &b_gm);
    qirp_inte_ligr(PT.delta, &a_dt, &b_dt);

    k = 0;
    step = (b_al - a_al) / (double) N;
    for (i = 0; i < N; i++) {
        t = a_al + (double) i *step;
        nehl_eval_segt(PT.alpha, t, &P[k]);
        k++;
    }

    step = (b_bt - a_bt) / (double) N;
    for (i = 0; i < N; i++) {
        t = a_bt + (double) i *step;
        nehl_eval_segt(PT.beta, t, &P[k]);
        k++;
    }

    step = (b_gm - a_gm) / (double) N;
    for (i = 0; i < N; i++) {
        t = a_gm + (double) i *step;
        nehl_eval_segt(PT.gamma, t, &P[k]);
        k++;
    }

    step = (b_dt - a_dt) / (double) N;
    for (i = 0; i < N; i++) {
        t = a_dt + (double) i *step;
        nehl_eval_segt(PT.delta, t, &P[k]);
        k++;
    }

    homs_boun_gosm(P, k, &B->x_min, &B->x_max, &B->y_min, &B->y_max, &B->z_min, &B->z_max);
    free(P);
    B->x_min = B->x_min - marg;
    B->x_max = B->x_max + marg;
    B->y_min = B->y_min - marg;
    B->y_max = B->y_max + marg;
    B->z_min = B->z_min - marg;
    B->z_max = B->z_max + marg;
}



void nojp_boun_zovq(pt_cnv PC, double marg, bd_box3D * B)
{
    int N = 5, i, k;
    double step, t;
    double a_al, b_al;
    double a_bt, b_bt;
    double a_gm, b_gm;
    point *P;

    P = (point *) malloc(3 * N * sizeof(point));
    qirp_inte_ligr(PC.alpha, &a_al, &b_al);
    qirp_inte_ligr(PC.beta, &a_bt, &b_bt);
    qirp_inte_ligr(PC.gamma, &a_gm, &b_gm);

    k = 0;
    step = (b_al - a_al) / (double) N;
    for (i = 0; i < N; i++) {
        t = a_al + (double) i *step;
        nehl_eval_segt(PC.alpha, t, &P[k]);
        k++;
    }

    step = (b_bt - a_bt) / (double) N;
    for (i = 0; i < N; i++) {
        t = a_bt + (double) i *step;
        nehl_eval_segt(PC.beta, t, &P[k]);
        k++;
    }

    step = (b_gm - a_gm) / (double) N;
    for (i = 0; i < N; i++) {
        t = a_gm + (double) i *step;
        nehl_eval_segt(PC.gamma, t, &P[k]);
        k++;
    }

    homs_boun_gosm(P, k, &B->x_min, &B->x_max, &B->y_min, &B->y_max, &B->z_min, &B->z_max);
    free(P);
    B->x_min = B->x_min - marg;
    B->x_max = B->x_max + marg;
    B->y_min = B->y_min - marg;
    B->y_max = B->y_max + marg;
    B->z_min = B->z_min - marg;
    B->z_max = B->z_max + marg;
}



double jufq_dist_dunt(side_flag SF, pt_tor PT, c_arc3D C, int *side, int *found_id)
{
    int i, q, *excl;
    double dis, res;
    c_arc3D *CA;
    CA = (c_arc3D *) malloc(4 * sizeof(c_arc3D));
    poms_find_resk_lonb(PT.alpha, &CA[0]);
    poms_find_resk_lonb(PT.beta, &CA[1]);
    poms_find_resk_lonb(PT.gamma, &CA[2]);
    poms_find_resk_lonb(PT.delta, &CA[3]);
    excl = (int *) malloc(4 * sizeof(int));
    if (SF.aflg == -1)
        excl[0] = 0;
    else
        excl[0] = 1;
    if (SF.bflg == -1)
        excl[1] = 0;
    else
        excl[1] = 1;
    if (SF.gflg == -1)
        excl[2] = 0;
    else
        excl[2] = 1;
    if (SF.dflg == -1)
        excl[3] = 0;
    else
        excl[3] = 1;
    res = LARGE_NUMBER;
    q = -1;
    *found_id = 0;
    for (i = 0; i < 4; i++)
        if (excl[i] == 0) {
            dis = fozm_dist_lojn(CA[i], C);
            if (dis < res) {
                res = dis;
                q = i;
                *found_id = 1;
            }
        }
    free(excl);
    free(CA);
    if (q == 0)
        *side = ON_ALPHA;
    if (q == 1)
        *side = ON_BETA;
    if (q == 2)
        *side = ON_GAMMA;
    if (q == 3)
        *side = ON_DELTA;
    return res;
}



double wujt_dist_kaqt(side_flag SF, pt_cnv PC, c_arc3D C, int *side, int *found_id)
{
    int i, q, *excl;
    double dis, res;
    c_arc3D *CA;
    CA = (c_arc3D *) malloc(3 * sizeof(c_arc3D));
    poms_find_resk_lonb(PC.alpha, &CA[0]);
    poms_find_resk_lonb(PC.beta, &CA[1]);
    poms_find_resk_lonb(PC.gamma, &CA[2]);
    excl = (int *) malloc(3 * sizeof(int));
    if (SF.aflg == -1)
        excl[0] = 0;
    else
        excl[0] = 1;
    if (SF.bflg == -1)
        excl[1] = 0;
    else
        excl[1] = 1;
    if (SF.gflg == -1)
        excl[2] = 0;
    else
        excl[2] = 1;
    res = LARGE_NUMBER;
    *found_id = 0;
    q = -1;
    for (i = 0; i < 3; i++)
        if (excl[i] == 0) {
            dis = fozm_dist_lojn(CA[i], C);
            if (dis < res) {
                res = dis;
                *found_id = 1;
                q = i;
            }
        }
    free(CA);
    free(excl);
    if (q == 0)
        *side = ON_ALPHA;
    if (q == 1)
        *side = ON_BETA;
    if (q == 2)
        *side = ON_GAMMA;
    return res;
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

