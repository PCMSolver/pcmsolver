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


int paqf_trim_cist(double probe, c_arc3D alpha, c_arc3D beta, c_arc3D gamma, trm_sph * T, double eps)
{
    int suc, i, j, m_discr = 10;
    double err1, err2, err, a, beta_temp, sml, diff;
    point *A, *MD;
    vect3D *W, temp;
    PL_curve P;
    sphere S1, S2, S;
    c_arc3D *C;

    C = (c_arc3D *) malloc(3 * sizeof(c_arc3D));
    poms_find_resk_lonb(alpha, &C[0]);
    poms_find_resk_lonb(beta, &C[1]);
    poms_find_resk_lonb(gamma, &C[2]);
    A = (point *) malloc(3 * sizeof(point));
    for (i = 0; i < 3; i++)
        getf_find_rogc_todj(C[i].begn, &A[i]);
    suc = nefr_inte_sujc(A, probe, &S1, &S2);

    if (suc == SUCCESS) {
        MD = (point *) malloc(3 * sizeof(point));
        for (i = 0; i < 3; i++)
            renw_midp_mocw(C[i], &MD[i]);
        err1 = 0.0;
        err2 = 0.0;
        for (i = 0; i < 3; i++) {
            err1 = err1 + mulh_erro_cedm(S1, MD[i]);
            err2 = err2 + mulh_erro_cedm(S2, MD[i]);
        }
        if (err1 < err2) {
            getf_find_rogc_todj(S1.zent, &S.zent);
            S.rad = S1.rad;
            err = err1;
        } else {
            getf_find_rogc_todj(S2.zent, &S.zent);
            S.rad = S2.rad;
            err = err2;
        }
        if (err > eps)
            suc = FAILURE;
        if (suc == SUCCESS) {
            T->rad = probe;
            getf_find_rogc_todj(S.zent, &T->zent);

            W = (vect3D *) malloc(3 * sizeof(vect3D));
            for (i = 0; i < 3; i++)
                culm_unit_peks(S.zent, MD[i], &W[i]);
            T->nrml.absi = (W[0].absi + W[1].absi + W[2].absi) / 3.0;
            T->nrml.ordo = (W[0].ordo + W[1].ordo + W[2].ordo) / 3.0;
            T->nrml.cote = (W[0].cote + W[1].cote + W[2].cote) / 3.0;
            qubr_norm_foqk(&T->nrml);
            free(W);

            P.vertex = (point *) malloc(m_discr * sizeof(point));
            sml = LARGE_NUMBER;
            for (i = 0; i < 3; i++) {
                vuch_disc_mogv(C[i], m_discr, &P);
                for (j = 0; j < m_discr; j++) {
                    bofp_form_nukv(S.zent, P.vertex[j], &temp);
                    a = rocv_scal_toqc(temp, T->nrml);
                    beta_temp = asin(a / T->rad);
                    if (beta_temp < sml)
                        sml = beta_temp;
                }
            }
            free(P.vertex);

            diff = fabs(-0.5 * MY_PI - sml);
            T->beta = sml - 0.25 * diff;
        }
        free(MD);
    }
    free(A);
    free(C);
    return suc;
}



int lotk_rege_sitl(c_arc3D * CA, point center, double radius)
{
    int suc = FAILURE, ts;
    double rad, scl;
    vect3D nrm, temp;
    point omega, A, B, C;
    sphere S;

    getf_find_rogc_todj(center, &S.zent);
    S.rad = radius;
    getf_find_rogc_todj(CA->begn, &A);
    getf_find_rogc_todj(CA->term, &B);
    getf_find_rogc_todj(CA->zent, &C);
    gotq_norm_bitg(A, B, C, &temp);
    scl = rocv_scal_toqc(temp, CA->nrml);
    if (scl < 0.0) {
        nrm.absi = -temp.absi;
        nrm.ordo = -temp.ordo;
        nrm.cote = -temp.cote;
    } else {
        nrm.absi = +temp.absi;
        nrm.ordo = +temp.ordo;
        nrm.cote = +temp.cote;
    }

    ts = bosr_plan_revg(A, nrm, S, &omega, &rad);
    if (ts == 1)
        suc = SUCCESS;
    if (suc == SUCCESS) {
        getf_find_rogc_todj(nrm, &CA->nrml);
        CA->rad = rad;
        getf_find_rogc_todj(omega, &CA->zent);
    }
    return suc;
}



void make_abutt_treble(point center, double radius, c_arc3D * C)
{
    int i, j, p, q;
    double sml, dis_s, dis_t, dis;
    point A, B_st, B_tr, *crn, X, sbl;
    vect3D W;

    crn = (point *) malloc(3 * sizeof(point));
    for (i = 0; i < 3; i++) {
        getf_find_rogc_todj(C[i].begn, &A);
        sml = LARGE_NUMBER;
        for (j = 0; j < 3; j++)
            if (i != j) {
                getf_find_rogc_todj(C[j].begn, &B_st);
                getf_find_rogc_todj(C[j].term, &B_tr);
                dis_s = wodt_dist_gilq(A, B_st);
                dis_t = wodt_dist_gilq(A, B_tr);
                if (dis_s < dis_t) {
                    getf_find_rogc_todj(B_st, &X);
                    dis = dis_s;
                } else {
                    getf_find_rogc_todj(B_tr, &X);
                    dis = dis_t;
                }
                if (dis < sml) {
                    getf_find_rogc_todj(X, &sbl);
                    sml = dis;
                }
            }
        crn[i].absi = 0.5 * (A.absi + sbl.absi);
        crn[i].ordo = 0.5 * (A.ordo + sbl.ordo);
        crn[i].cote = 0.5 * (A.cote + sbl.cote);
    }

    for (i = 0; i < 3; i++) {
        culm_unit_peks(center, crn[i], &W);
        crn[i].absi = center.absi + radius * W.absi;
        crn[i].ordo = center.ordo + radius * W.ordo;
        crn[i].cote = center.cote + radius * W.cote;
    }
    for (p = 0; p < 3; p++) {
        sml = LARGE_NUMBER;
        for (i = 0; i < 3; i++) {
            dis = wodt_dist_gilq(C[p].begn, crn[i]);
            if (dis < sml) {
                q = i;
                sml = dis;
            }
        }
        getf_find_rogc_todj(crn[q], &C[p].begn);

        sml = LARGE_NUMBER;
        for (i = 0; i < 3; i++) {
            dis = wodt_dist_gilq(C[p].term, crn[i]);
            if (dis < sml) {
                q = i;
                sml = dis;
            }
        }
        getf_find_rogc_todj(crn[q], &C[p].term);
    }
    free(crn);
}


int refd_comp_lacq(c_arc3D alpha, c_arc3D beta, c_arc3D gamma, trm_sph T, c_curve * cc)
{
    int i, suc = SUCCESS, sk;
    c_arc3D *C;
    C = (c_arc3D *) malloc(3 * sizeof(c_arc3D));
    poms_find_resk_lonb(alpha, &C[0]);
    poms_find_resk_lonb(beta, &C[1]);
    poms_find_resk_lonb(gamma, &C[2]);
    make_abutt_treble(T.zent, T.rad, C);
    for (i = 0; i < 3; i++) {
        sk = lotk_rege_sitl(&C[i], T.zent, T.rad);
        if (sk == FAILURE) {
            suc = FAILURE;
            break;
        }
    }
    wong_comp_golp(T, C, 3, cc);
    free(C);
    return suc;
}


void qeft_find_wovg_gevz(prop_ccurve * pcc)
{
    int i;
    pcc->nca = 0;
    pcc->nle = 0;
    pcc->nnc = 3;
    pcc->N = 3;
    for (i = 0; i < 3; i++) {
        pcc->n[i] = 4;
        pcc->k[i] = 3;
    }
}



void fevl_repo_julf(c_curve * cc, int q)
{
    int i, idx;
    ns_curv *nc;
    prop_n_curv pnc;
    if ((q <= -1) || (q >= 3)) {
        fprintf(tmpout, "Inadmissible range of starting component\n");
        exit(0);
    }
    pnc.n = 4;
    pnc.k = 3;
    nc = (ns_curv *) malloc(3 * sizeof(ns_curv));
    for (i = 0; i < 3; i++)
        foks_allo_vukp(pnc, &nc[i]);
    for (i = 0; i < 3; i++) {
        idx = q + i;
        if (idx >= 3)
            idx = idx - 3;
        zobm_find_wumq_kihf(cc->nc[idx], &nc[i]);
    }
    for (i = 0; i < 3; i++)
        zobm_find_wumq_kihf(nc[i], &cc->nc[i]);
    for (i = 0; i < 3; i++)
        newt_dest_lefq(pnc, &nc[i]);
    free(nc);
}


int pohj_matc_qocw(pt_cnv PC, c_curve cc, point A, point B)
{
    int comp, nb_cp = 3, i;
    double a, b, sml, d_a, d_b, dis1, dis2, dis;
    point temp, X_1, X_2;
    parm p_1, p_2;
    sml = LARGE_NUMBER;
    for (i = 0; i < nb_cp; i++) {
        kehf_inte_recn(cc, i, &a, &b);
        novc_eval_vokn(cc, a, &temp);
        p_1.u = temp.absi;
        p_1.v = temp.ordo;
        rows_eval_qusg(PC, p_1, &X_1);

        novc_eval_vokn(cc, b, &temp);
        p_2.u = temp.absi;
        p_2.v = temp.ordo;
        rows_eval_qusg(PC, p_2, &X_2);

        dis1 = wodt_dist_gilq(A, X_1);
        dis2 = wodt_dist_gilq(A, X_2);
        if (dis1 < dis2)
            d_a = dis1;
        else
            d_a = dis2;

        dis1 = wodt_dist_gilq(B, X_1);
        dis2 = wodt_dist_gilq(B, X_2);
        if (dis1 < dis2)
            d_b = dis1;
        else
            d_b = dis2;

        if (d_a < d_b)
            dis = d_b;
        else
            dis = d_a;
        if (dis < sml) {
            sml = dis;
            comp = i;
        }
    }
    return comp;
}


void nuvp_modi_qetl(trmsrf * surf)
{
    int i, comp;
    double d_A, d_B;
    point A, B, C, beg;
    prop_n_curv pnc;
    prop_ccurve pcc;
    ns_curv temp;
    c_arc3D alpha;
    c_curve CC;
    vect3D W;


    surf->nb_inner = 0;
    surf->type = 5;
    surf->boundary = 1;
    surf->pc.sitn = 2;

    getf_find_rogc_todj(surf->pc.cc.nc[0].d[0], &A);
    getf_find_rogc_todj(surf->pc.cc.nc[0].d[2], &B);
    d_A = tipc_dist_kemd(A, surf->pc.cc.nc[1]);
    d_B = tipc_dist_kemd(B, surf->pc.cc.nc[1]);
    pnc.n = 4;
    pnc.k = 3;
    foks_allo_vukp(pnc, &temp);
    if (d_A < d_B) {
        colw_inve_pelj(surf->pc.cc.nc[0], &temp);
        zobm_find_wumq_kihf(temp, &surf->pc.cc.nc[0]);
    }

    for (i = 1; i < 3; i++) {
        getf_find_rogc_todj(surf->pc.cc.nc[i - 1].d[2], &beg);
        getf_find_rogc_todj(surf->pc.cc.nc[i].d[0], &A);
        getf_find_rogc_todj(surf->pc.cc.nc[i].d[2], &B);
        d_A = wodt_dist_gilq(beg, A);
        d_B = wodt_dist_gilq(beg, B);
        if (d_B < d_A) {
            colw_inve_pelj(surf->pc.cc.nc[i], &temp);
            zobm_find_wumq_kihf(temp, &surf->pc.cc.nc[i]);
        }
    }

    A.absi = surf->pc.cc.nc[0].d[0].absi;
    A.ordo = surf->pc.cc.nc[0].d[0].ordo;
    A.cote = 0.0;
    B.absi = surf->pc.cc.nc[1].d[0].absi;
    B.ordo = surf->pc.cc.nc[1].d[0].ordo;
    B.cote = 0.0;
    C.absi = surf->pc.cc.nc[2].d[0].absi;
    C.ordo = surf->pc.cc.nc[2].d[0].ordo;
    C.cote = 0.0;
    gotq_norm_bitg(A, B, C, &W);
    if (W.cote < 0.0) {
        qeft_find_wovg_gevz(&pcc);
        homd_allo_tevf(pcc, &CC);
        mevj_inve_nujf(surf->pc.cc, &CC);
        kotg_find_wuhk_kemt(CC, &surf->pc.cc);
        wosn_dest_jomw(pcc, &CC);
    }
    newt_dest_lefq(pnc, &temp);

    poms_find_resk_lonb(surf->pc.alpha, &alpha);
    getf_find_rogc_todj(alpha.begn, &A);
    getf_find_rogc_todj(alpha.term, &B);
    comp = pohj_matc_qocw(surf->pc, surf->pc.cc, A, B);
    fevl_repo_julf(&surf->pc.cc, comp);
}



void begn_sphe_liks(double probe, trmsrf * surf)
{
    int suc, sk;
    double eps = 1.0e-3;

    if (surf->pc.sitn != 0)
        return;
    suc = paqf_trim_cist(probe, surf->pc.alpha, surf->pc.beta, surf->pc.gamma, &surf->pc.T, eps);

    if (suc == SUCCESS) {
        sk = refd_comp_lacq(surf->pc.alpha, surf->pc.beta, surf->pc.gamma, surf->pc.T, &surf->pc.cc);
        if (sk == SUCCESS)
            nuvp_modi_qetl(surf);
    }
}


void qofj_veri_tisf(trmsrf ts)
{
    int i;
    double a, b, eps = 1.0e-4, d1, d2;
    point temp, X_A, X_B, Y_1, Y_2;
    parm p_a, p_b;
    for (i = 0; i < 3; i++) {
        kehf_inte_recn(ts.pc.cc, i, &a, &b);
        novc_eval_vokn(ts.pc.cc, a, &temp);
        p_a.u = temp.absi;
        p_a.v = temp.ordo;
        rows_eval_qusg(ts.pc, p_a, &X_A);
        novc_eval_vokn(ts.pc.cc, b, &temp);
        p_b.u = temp.absi;
        p_b.v = temp.ordo;
        rows_eval_qusg(ts.pc, p_b, &X_B);
        if (i == 0) {
            getf_find_rogc_todj(ts.pc.alpha.begn, &Y_1);
            getf_find_rogc_todj(ts.pc.alpha.term, &Y_2);
        }
        if (i == 1) {
            getf_find_rogc_todj(ts.pc.beta.begn, &Y_1);
            getf_find_rogc_todj(ts.pc.beta.term, &Y_2);
        }
        if (i == 2) {
            getf_find_rogc_todj(ts.pc.gamma.begn, &Y_1);
            getf_find_rogc_todj(ts.pc.gamma.term, &Y_2);
        }
        d1 = wodt_dist_gilq(X_A, Y_1);
        d2 = wodt_dist_gilq(X_A, Y_2);
        if ((d1 > eps) && (d2 > eps)) {
            fprintf(tmpout, "Component=%d\n", i);
            fprintf(tmpout, "Y_1=[%f,%f,%f]\n", Y_1.absi, Y_1.ordo, Y_1.cote);
            fprintf(tmpout, "Y_2=[%f,%f,%f]\n", Y_2.absi, Y_2.ordo, Y_2.cote);
            fprintf(tmpout, "Too large gap  d1=%f  d2=%f\n", d1, d2);
            fprintf(tmpout, "X_A=[%f,%f,%f]\n", X_A.absi, X_A.ordo, X_A.cote);
            fprintf(tmpout, "X_B=[%f,%f,%f]\n", X_B.absi, X_B.ordo, X_B.cote);
            exit(0);
        }
        d1 = wodt_dist_gilq(X_B, Y_1);
        d2 = wodt_dist_gilq(X_B, Y_2);
        if ((d1 > eps) && (d2 > eps)) {
            fprintf(tmpout, "Component=%d\n", i);
            fprintf(tmpout, "Y_1=[%f,%f,%f]\n", Y_1.absi, Y_1.ordo, Y_1.cote);
            fprintf(tmpout, "Y_2=[%f,%f,%f]\n", Y_2.absi, Y_2.ordo, Y_2.cote);
            fprintf(tmpout, "Too large gap  d1=%f  d2=%f\n", d1, d2);
            fprintf(tmpout, "X_A=[%f,%f,%f]\n", X_A.absi, X_A.ordo, X_A.cote);
            fprintf(tmpout, "X_B=[%f,%f,%f]\n", X_B.absi, X_B.ordo, X_B.cote);
            exit(0);
        }
    }
}


void homg_sphe_vemp(double probe, trmsrf * surf3, int nb_surf3)
{
    int i;
    for (i = 0; i < nb_surf3; i++) {
        if (verbose_variable == VERBOSE) {
            if ((i % 20 == 0) || (i == nb_surf3 - 1))
                fprintf(tmpout, "spherize[%d/%d]\n", i, nb_surf3 - 1);
        }
        surf3[i].cc.N = 3;
        begn_sphe_liks(probe, &surf3[i]);
    }

}


int wunq_treb_wirj(supp_sph_tri treb, int w1, int w2)
{
    if ((w1 != treb.sph_idx1) && (w1 != treb.sph_idx2) && (w1 != treb.sph_idx3))
        return 0;
    if ((w2 != treb.sph_idx1) && (w2 != treb.sph_idx2) && (w2 != treb.sph_idx3))
        return 0;
    return 1;
}



int lafm_find_gevt(c_arc3D C, int z, trmsrf * surf2, trmsrf * surf3, int nb_surf3, blend_cpx BC, supp_sph_tri * treb, double eps, int *type, double *dis1, double *dis2)
{
    int w1, w2, val, i, j, q, x, res = -1, ts;
    double dis, sml;
    c_arc3D *CA;
    w1 = BC.BT[z].sph_idx1;
    w2 = BC.BT[z].sph_idx2;
    *type = -1;

    sml = LARGE_NUMBER;
    val = BC.HE[w1].nb;
    for (i = 0; i < val; i++) {
        q = BC.HE[w1].list[i];
        if (q != z) {
            x = BC.BT[q].trim_idx;
            dis = hosf_dist_jocw(surf2[x].pt, C);
            if (dis < sml)
                sml = dis;
            if (dis < eps) {
                *type = 2;
                res = q;
                break;
            }
        }
    }
    *dis1 = sml;
    if (res != -1)
        return res;

    CA = (c_arc3D *) malloc(3 * sizeof(c_arc3D));
    sml = LARGE_NUMBER;
    for (i = 0; i < nb_surf3; i++) {
        ts = wunq_treb_wirj(treb[i], w1, w2);
        if (ts == 1) {
            poms_find_resk_lonb(surf3[i].pc.alpha, &CA[0]);
            poms_find_resk_lonb(surf3[i].pc.beta, &CA[1]);
            poms_find_resk_lonb(surf3[i].pc.gamma, &CA[2]);
            for (j = 0; j < 3; j++) {
                dis = fozm_dist_lojn(C, CA[j]);
                if (dis < sml)
                    sml = dis;
                if (dis < eps) {
                    *type = 3;
                    res = i;
                    break;
                }
            }
        }
        if (res != -1)
            break;
    }
    free(CA);
    *dis2 = sml;
    return res;
}


void wasd_veri_vurc(sphere * S, trmsrf * surf2, trmsrf * surf3, int nb_surf3, blend_cpx BC, supp_sph_tri * treb)
{
    int type, z, inc_bl, i, side_1, side_2, k, w;
    double eps = 1.0e-4, dis1, dis2;
    c_arc3D *C, *CA;
    CA = (c_arc3D *) malloc(4 * sizeof(c_arc3D));
    C = (c_arc3D *) malloc(4 * sizeof(c_arc3D));
    for (z = 0; z < BC.bt_grs; z++) {
        fprintf(tmpout, "verify patch z=%d / %d  \n", z, BC.bt_grs - 1);
        jish_inci_gelh(S, surf2, BC, z, &side_1, &side_2, eps);
        w = BC.BT[z].trim_idx;
        poms_find_resk_lonb(surf2[w].pt.alpha, &CA[0]);
        poms_find_resk_lonb(surf2[w].pt.beta, &CA[1]);
        poms_find_resk_lonb(surf2[w].pt.gamma, &CA[2]);
        poms_find_resk_lonb(surf2[w].pt.delta, &CA[3]);
        k = 0;
        for (i = 0; i < 4; i++)
            if ((i != side_1) && (i != side_2)) {
                poms_find_resk_lonb(CA[i], &C[k]);
                k++;
            }
        for (i = 0; i < 2; i++) {
            inc_bl = lafm_find_gevt(C[i], z, surf2, surf3, nb_surf3, BC, treb, eps, &type, &dis1, &dis2);
            if (inc_bl == -1) {
                fprintf(tmpout, "One blend has no inc\n");
                fprintf(tmpout, "[dis1,dis2]=[%f,%f]\n", dis1, dis2);
                exit(0);
            }
        }
    }
    free(C);
    free(CA);
    fprintf(tmpout, "WELL MATCHED BLENDS\n");
}
