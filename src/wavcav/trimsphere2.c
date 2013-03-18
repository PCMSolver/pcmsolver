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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"
#include "sas.h"


void qatf_conv_ducq(ns_curv nc, c_arc2D * CA)
{
    int n, i, N = 20;
    double rad, a, b, tm, dis0, dis1, D0, D1;
    parm omega, S, M, T, *disc0, *disc1;
    c_arc2D *cand;
    S.u = nc.d[0].absi;
    S.v = nc.d[0].ordo;
    n = nc.n;
    T.u = nc.d[n].absi;
    T.v = nc.d[n].ordo;
    a = nc.v0;
    b = nc.v1;
    tm = 0.5 * (a + b);
    kuqt_eval_webp(nc, tm, &M);
    cehk_circ_jesw(S, M, T, &rad, &omega);

    cand = (c_arc2D *) malloc(2 * sizeof(c_arc2D));
    for (i = 0; i < 2; i++) {
        cunl_find_qedf_rewn(omega, &cand[i].zent);
        cand[i].c_cir = 0;
        cand[i].rad = rad;
    }
    cunl_find_qedf_rewn(S, &cand[0].begn);
    cunl_find_qedf_rewn(T, &cand[0].term);
    cunl_find_qedf_rewn(T, &cand[1].begn);
    cunl_find_qedf_rewn(S, &cand[1].term);

    disc0 = (parm *) malloc(N * sizeof(parm));
    disc1 = (parm *) malloc(N * sizeof(parm));
    lohs_disc_gojz(cand[0], N, disc0);
    lohs_disc_gojz(cand[1], N, disc1);
    dis0 = LARGE_NUMBER;
    dis1 = LARGE_NUMBER;
    for (i = 0; i < N; i++) {
        D0 = pufv_dist_mekq(disc0[i], M);
        D1 = pufv_dist_mekq(disc1[i], M);
        if (D0 < dis0)
            dis0 = D0;
        if (D1 < dis1)
            dis1 = D1;
    }
    free(disc0);
    free(disc1);

    if (dis0 < dis1)
        goth_find_cofl_futw(cand[0], CA);
    else
        goth_find_cofl_futw(cand[1], CA);
    free(cand);
}


void kugw_conv_nimq(ns_curv C, c_arc2D * CA)
{
    int N = 10, i;
    double xmi, xma, ymi, yma, hx, hy, lrg;
    c_arc2D CA_aux;
    ns_curv C_aux;
    prop_n_curv pnc;
    parm *p, G;

    p = (parm *) malloc(N * sizeof(parm));
    tegn_disc_likp(C, N, p);
    ritp_boun_niwz(p, N, &xmi, &xma, &ymi, &yma);
    G.u = 0.5 * (xmi + xma);
    G.v = 0.5 * (ymi + yma);
    hx = xma - xmi;
    hy = yma - ymi;
    if (hx > hy)
        lrg = hx;
    else
        lrg = hy;
    free(p);

    pnc.k = C.k;
    pnc.n = C.n;
    foks_allo_vukp(pnc, &C_aux);
    zobm_find_wumq_kihf(C, &C_aux);
    for (i = 0; i <= C.n; i++) {
        C_aux.d[i].absi = (C.d[i].absi - G.u) / lrg;
        C_aux.d[i].ordo = (C.d[i].ordo - G.v) / lrg;
        C_aux.d[i].cote = 0.0;
    }
    qatf_conv_ducq(C_aux, &CA_aux);
    newt_dest_lefq(pnc, &C_aux);

    goth_find_cofl_futw(CA_aux, CA);
    CA->zent.u = (lrg * CA_aux.zent.u) + G.u;
    CA->zent.v = (lrg * CA_aux.zent.v) + G.v;
    CA->begn.u = (lrg * CA_aux.begn.u) + G.u;
    CA->begn.v = (lrg * CA_aux.begn.v) + G.v;
    CA->term.u = (lrg * CA_aux.term.u) + G.u;
    CA->term.v = (lrg * CA_aux.term.v) + G.v;
    CA->rad = lrg * CA_aux.rad;
}



int lepg_find_rijf_fenv(c_arc2D C, double x, double y, parm * I)
{
    int i, ind, suc = FAILURE;
    double temp, theta, alpha_s, alpha_t;
    double *cand, alpha;
    temp = y / C.rad;
    theta = asin(temp);
    cand = (double *) malloc(2 * sizeof(double));
    cand[0] = +theta;
    cand[1] = MY_PI - theta;
    laqv_inte_jotl(C, &alpha_s, &alpha_t);
    for (i = 0; i < 2; i++) {
        alpha = cand[i];
        ind = 0;
        if ((alpha_s < alpha) && (alpha < alpha_t))
            ind = 1;
        if (ind == 0)
            if ((alpha_s < alpha + 2.0 * MY_PI) && (alpha + 2.0 * MY_PI < alpha_t))
                ind = 1;
        if (ind == 1) {
            I->u = C.rad * cos(alpha);
            I->v = C.rad * sin(alpha);
            suc = SUCCESS;
            break;
        }
    }
    free(cand);
    return suc;
}



void ponv_stat_jozq(c_arc2D C, int *state_bt, int *state_tp, double *al_bt, double *al_tp)
{
    int i;
    double a, b, demi_pi, trs_pi;
    laqv_inte_jotl(C, &a, &b);
    *state_tp = 0;
    demi_pi = 0.5 * MY_PI;
    for (i = 0; i < 3; i++) {
        if ((a < demi_pi) && (demi_pi < b)) {
            *state_tp = 1;
            *al_tp = demi_pi;
            break;
        }
        demi_pi = demi_pi + 2.0 * MY_PI;
    }

    *state_bt = 0;
    trs_pi = -0.5 * MY_PI;
    for (i = 0; i < 3; i++) {
        if ((a < trs_pi) && (trs_pi < b)) {
            *state_bt = 1;
            *al_bt = trs_pi;
            break;
        }
        trs_pi = trs_pi + 2.0 * MY_PI;
    }
}


void zacd_spli_tugd(c_arc2D C, double *al, int nb, c_arc2D * c)
{
    int i;
    parm *sep;
    if (nb == 0) {
        goth_find_cofl_futw(C, &c[0]);
        return;
    }
    sep = (parm *) malloc((nb + 2) * sizeof(parm));
    cunl_find_qedf_rewn(C.begn, &sep[0]);
    for (i = 0; i < nb; i++)
        logm_eval_pusn(C, al[i], &sep[i + 1]);
    cunl_find_qedf_rewn(C.term, &sep[nb + 1]);
    for (i = 0; i <= nb; i++) {
        goth_find_cofl_futw(C, &c[i]);
        cunl_find_qedf_rewn(sep[i], &c[i].begn);
        cunl_find_qedf_rewn(sep[i + 1], &c[i].term);
    }
    free(sep);
}


int mopj_spli_lerw(c_arc2D C, c_arc2D * c)
{
    int st_bt, st_tp, nb, nb_sep;
    double al_bt, al_tp, *al, temp;
    ponv_stat_jozq(C, &st_bt, &st_tp, &al_bt, &al_tp);
    if ((st_bt == 0) && (st_tp == 0)) {
        goth_find_cofl_futw(C, &c[0]);
        nb = 1;
    } else {
        al = (double *) malloc(2 * sizeof(double));
        nb_sep = 0;
        if (st_tp == 1) {
            al[nb_sep] = al_tp;
            nb_sep++;
        }
        if (st_bt == 1) {
            al[nb_sep] = al_bt;
            nb_sep++;
        }

        if (nb_sep == 2) {
            if (al[0] > al[1]) {
                temp = al[0];
                al[0] = al[1];
                al[1] = temp;
            }
        }

        zacd_spli_tugd(C, al, nb_sep, c);
        nb = nb_sep + 1;
        free(al);
    }
    return nb;
}


int vetc_prep_fems(c_arc2D * C_in, int N, c_arc2D * C_out)
{
    int i, j, M, nb;
    c_arc2D *temp;
    temp = (c_arc2D *) malloc(4 * sizeof(c_arc2D));
    M = 0;
    for (i = 0; i < N; i++) {
        nb = mopj_spli_lerw(C_in[i], temp);
        for (j = 0; j < nb; j++) {
            goth_find_cofl_futw(temp[j], &C_out[M]);
            M++;
        }
    }
    free(temp);
    return M;
}


void hegf_span_taqh(c_arc2D C, double *y_mn, double *y_mx)
{
    double cand_st, cand_tr, sml, lrg;
    cand_st = C.begn.v;
    cand_tr = C.term.v;
    if (cand_st < cand_tr) {
        sml = cand_st;
        lrg = cand_tr;
    } else {
        sml = cand_tr;
        lrg = cand_st;
    }
    *y_mn = sml;
    *y_mx = lrg;
}



int gulj_arco_hivf(double y, double eps, c_arc2D C)
{
    double y_d, y_u, diff, y_mn, y_mx;

    y_d = C.zent.v - C.rad;
    if (y < y_d)
        return 0;

    y_u = C.zent.v + C.rad;
    if (y > y_u)
        return 0;

    diff = fabs(y - C.begn.v);
    if (diff < eps)
        return 2;

    diff = fabs(y - C.term.v);
    if (diff < eps)
        return 2;

    hegf_span_taqh(C, &y_mn, &y_mx);
    if ((y_mn < y) && (y < y_mx))
        return 1;
    return 0;
}



int wilz_arco_lokn(double x, double y, double eps, c_arc2D C)
{
    int stn, suc;
    double x_aux, y_aux;
    parm I;
    c_arc2D C_aux;
    stn = gulj_arco_hivf(y, eps, C);

    if ((stn == 0) || (stn == 2))
        return stn;
    goth_find_cofl_futw(C, &C_aux);
    C_aux.zent.u = 0.0;
    C_aux.zent.v = 0.0;
    C_aux.begn.u = C.begn.u - C.zent.u;
    C_aux.begn.v = C.begn.v - C.zent.v;
    C_aux.term.u = C.term.u - C.zent.u;
    C_aux.term.v = C.term.v - C.zent.v;
    x_aux = x - C.zent.u;
    y_aux = y - C.zent.v;
    suc = lepg_find_rijf_fenv(C_aux, x_aux, y_aux, &I);
    if (suc == FAILURE)
        return 0;

    if (I.u < x_aux)
        return 0;
    return 1;
}



int vahg_find_purt_tumh(parm p, c_arc2D * C, int N)
{
    int ts, i, res, nb, par;
    double eps = 1.0e-4;
    nb = 0;
    for (i = 0; i < N; i++) {
        ts = wilz_arco_lokn(p.u, p.v, eps, C[i]);

        if (ts == 2) {
            res = INTERF_SET;
            break;
        }
        if (ts == 1)
            nb++;
    }

    par = nb % 2;
    if (par == 1)
        res = INS_SET;
    if (par == 0)
        res = OUT_SET;
    return res;
}


void wudg_find_zoml_dafw(polyarc * cc, int n, int N, mult_conn * mc)
{
    int i, j, k, l, m;
    parm *p;
    mc->nb_inner_polygons = n - 1;
    mc->nb_vr_outer = N * cc[0].ar_grs;
    for (i = 1; i < n; i++)
        mc->nb_vr_inner[i - 1] = N * cc[i].ar_grs;
    p = (parm *) malloc(N * sizeof(parm));
    k = 0;
    for (i = 0; i < n; i++) {
        m = cc[i].ar_grs;
        for (j = 0; j < m; j++) {
            lohs_disc_gojz(cc[i].CA[j], N, p);
            for (l = 0; l < N; l++) {
                cunl_find_qedf_rewn(p[l], &mc->vertex[k]);
                k++;
            }
        }
    }
    free(p);
    mc->v_grs = k;
}


void hogl_find_jipq_magv(polyarc pa_in, polyarc * pa_out)
{
    int N, i;
    N = pa_in.ar_grs;
    for (i = 0; i < N; i++)
        goth_find_cofl_futw(pa_in.CA[i], &pa_out->CA[i]);
    pa_out->ar_grs = N;
}


void lohs_disc_gojz(c_arc2D C, int N, parm * p)
{
    int i;
    double a, b, step, t;
    parm temp;
    laqv_inte_jotl(C, &a, &b);
    step = (b - a) / ((double) N - 1.0);
    for (i = 0; i < N; i++) {
        t = a + (double) i *step;
        logm_eval_pusn(C, t, &temp);
        p[i].u = temp.u;
        p[i].v = temp.v;
    }
}


void comf_bbox_renb(c_arc2D C, bd_box2D * BX)
{
    int N = 5;
    double xmi, xma, ymi, yma;
    parm *p;
    p = (parm *) malloc(N * sizeof(parm));
    lohs_disc_gojz(C, 3, p);
    vusf_clou_kegr(p, N, &xmi, &xma, &ymi, &yma);
    free(p);
    BX->x_min = xmi;
    BX->x_max = xma;
    BX->y_min = ymi;
    BX->y_max = yma;
}


void cegd_bbox_cosm(polyarc A, bd_box2D * BX)
{
    int i;
    double xmi, xma, ymi, yma;
    bd_box2D bx;
    xmi = LARGE_NUMBER;
    xma = -LARGE_NUMBER;
    ymi = LARGE_NUMBER;
    yma = -LARGE_NUMBER;
    for (i = 0; i < A.ar_grs; i++) {
        comf_bbox_renb(A.CA[i], &bx);
        if (bx.x_min < xmi)
            xmi = bx.x_min;
        if (bx.x_max > xma)
            xma = bx.x_max;
        if (bx.y_min < ymi)
            ymi = bx.y_min;
        if (bx.y_max > yma)
            yma = bx.y_max;
    }
    BX->x_min = xmi;
    BX->x_max = xma;
    BX->y_min = ymi;
    BX->y_max = yma;
}





int qiwg_find_johq_logk(polyarc W, bd_box2D BX)
{
    int res, i, j;
    parm *P;
    P = (parm *) malloc(2 * sizeof(parm));
    res = 1;
    for (i = 0; i < W.ar_grs; i++) {
        cunl_find_qedf_rewn(W.CA[i].begn, &P[0]);
        cunl_find_qedf_rewn(W.CA[i].term, &P[1]);
        for (j = 0; j < 2; j++) {
            if ((P[j].u < BX.x_min) || (P[j].u > BX.x_max)) {
                res = 0;
                break;
            }
            if ((P[j].v < BX.y_min) || (P[j].v > BX.y_max)) {
                res = 0;
                break;
            }
        }
        if (res == 0)
            break;
    }
    free(P);
    return res;
}



int mupv_chec_nefg(polyarc A, polyarc B)
{
    int res, N_A, N_B, pos, i, j, m_loc = 4, ts_bx;
    double a, b, t, step, lambda;
    parm p;
    bd_box2D BX;
    cegd_bbox_cosm(B, &BX);
    ts_bx = qiwg_find_johq_logk(A, BX);
    if (ts_bx == 0)
        return 0;
    N_A = A.ar_grs;
    N_B = B.ar_grs;
    step = 1.0 / (double) m_loc;
    res = -1;
    for (i = 0; i < N_A; i++) {
        laqv_inte_jotl(A.CA[i], &a, &b);
        for (j = 0; j < m_loc; j++) {
            lambda = (double) j *step;
            t = (1.0 - lambda) * a + lambda * b;
            logm_eval_pusn(A.CA[i], t, &p);
            pos = vahg_find_purt_tumh(p, B.CA, N_B);
            if (pos != INTERF_SET) {
                if (pos == INS_SET)
                    res = 1;
                if (pos == OUT_SET)
                    res = 0;
                break;
            }
        }
        if ((res == 0) || (res == 1))
            break;
    }

    return res;
}


void hotc_disp_tifl(polyarc P)
{
    int N, i;
    N = P.ar_grs;
    fprintf(tmpout, "Number of arcs=%d\n", N);
    for (i = 0; i < N; i++) {
        lutc_disp_nulq(P.CA[i]);
        fprintf(tmpout, "---\n");
    }
}



int mosc_encl_hapk(polyarc * P, int N, int k, int *excl, int m, int *list)
{
    int nb, i, ts, tr, dummy;
    nb = 0;
    for (i = 0; i < N; i++)
        if (i != k) {
            ts = gonl_arra_govj(excl, m, i, &dummy);
            if (ts == 0) {
                tr = mupv_chec_nefg(P[i], P[k]);
                if (tr == 1) {
                    list[nb] = i;
                    nb++;
                }
            }
        }
    return nb;
}



void pomq_curr_jecf(polygon * P, polyarc * P_arc, int n, int *excl, int m, int *EX, int *IN, int *NIN)
{
    int ex, nin, *list, i;
    ex = tulj_curr_winq(P, n, excl, m);
    list = (int *) malloc(n * sizeof(int));
    if (ex != -1) {
        nin = mosc_encl_hapk(P_arc, n, ex, excl, m, list);
        *NIN = nin;
        *EX = ex;
        for (i = 0; i < nin; i++)
            IN[i] = list[i];
    }
    free(list);
}


void jets_find_cifl_kitq(c_curve C, polyarc * P)
{
    int i, j, k, M, m;
    c_arc2D *CA, *CA_new;
    m = C.N;
    CA = (c_arc2D *) malloc(m * sizeof(c_arc2D));
    for (j = 0; j < m; j++) {
        kugw_conv_nimq(C.nc[j], &CA[j]);

    }
    CA_new = (c_arc2D *) malloc(4 * m * sizeof(c_arc2D));
    M = vetc_prep_fems(CA, m, CA_new);
    free(CA);
    k = 0;
    for (i = 0; i < M; i++) {
        goth_find_cofl_futw(CA_new[i], &P->CA[k]);
        k++;
    }
    P->ar_grs = k;
    free(CA_new);


}



int qonk_form_sevt(int par, trm_sph T, int nb_arcs, c_curve * cc, int n, trmsrf * surf, int nb_cur, int *supp, int max_surf)
{
    int NB, N = 20, i, m, *excl, nin, ex, *in, s, *orient;
    polygon *Q;
    polyarc *Q_arc;

    Q = (polygon *) malloc(n * sizeof(polygon));
    Q_arc = (polyarc *) malloc(n * sizeof(polyarc));
    for (i = 0; i < n; i++) {
        m = cc[i].N;
        Q[i].vertex = (parm *) malloc(m * N * sizeof(parm));
        Q[i].nb_local_vertices = (int *) malloc(sizeof(int));
        Q_arc[i].CA = (c_arc2D *) malloc(4 * m * sizeof(c_arc2D));
        jets_find_cifl_kitq(cc[i], &Q_arc[i]);
        fepm_find_jemr_wozg(cc[i], N, &Q[i]);


    }

    excl = (int *) malloc(n * sizeof(int));
    in = (int *) malloc(n * sizeof(int));
    m = 0;
    NB = 0;

    for (s = 0; s < n; s++) {
        pomq_curr_jecf(Q, Q_arc, n, excl, m, &ex, in, &nin);
        if (ex == -1) {
            fprintf(tmpout, "Unable to find exterior composite curve\n");
            exit(0);
        }
        orient = (int *) malloc((n + 1) * sizeof(int));
        if (nb_cur + NB >= max_surf) {
            fprintf(tmpout, "maximum number of trimmed surfaces is reached\n");
            exit(0);
        }
        pugj_crea_sotm(ex, in, nin, T, cc, n, &surf[nb_cur + NB], orient);
        supp[nb_cur + NB] = par;
        NB++;
        free(orient);

        excl[m] = ex;
        for (i = 0; i < nin; i++) {
            if (m + i + 1 > n) {
                fprintf(tmpout, "index exceeded\n");
                exit(0);
            }
            excl[m + i + 1] = in[i];
        }
        m = m + nin + 1;
        if (m >= n)
            break;
    }

    free(in);
    free(excl);
    for (i = 0; i < n; i++) {
        free(Q[i].vertex);
        free(Q[i].nb_local_vertices);
        free(Q_arc[i].CA);
    }
    free(Q);
    free(Q_arc);
    return NB;
}
