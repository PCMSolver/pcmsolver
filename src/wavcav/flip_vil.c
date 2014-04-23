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
#include "geodesic.h"
#include "pln_sph.h"
#include "smooth.h"


int nigw_peri_neql(int q)
{
    int res;
    res = q + 1;
    if (res >= 6)
        res = res - 6;
    return res;
}


void bugq_assi_sokm(int val, int w, int id, fajor_sion3D * quad)
{
    if (id == 1)
        quad->elem[w].frvrt = val;
    if (id == 2)
        quad->elem[w].scvrt = val;
    if (id == 3)
        quad->elem[w].thvrt = val;
    if (id == 4)
        quad->elem[w].ftvrt = val;
}



int vogj_edge_tocf(int n_a, int n_b, fajor_sion3D quad, int E, int F, int *id)
{
    int *cand, n1, n2, i, res = -1, e;
    cand = (int *) malloc(8 * sizeof(int));
    cand[0] = quad.elem[E].frkt;
    cand[1] = quad.elem[E].sckt;
    cand[2] = quad.elem[E].trkt;
    cand[3] = quad.elem[E].ftkt;
    cand[4] = quad.elem[F].frkt;
    cand[5] = quad.elem[F].sckt;
    cand[6] = quad.elem[F].trkt;
    cand[7] = quad.elem[F].ftkt;
    for (i = 0; i < 8; i++) {
        e = cand[i];
        n1 = quad.kt[e].frvrt;
        n2 = quad.kt[e].scvrt;
        if (((n1 == n_a) && (n2 == n_b)) || ((n2 == n_a) && (n1 == n_b))) {
            res = e;
            if (i < 4)
                *id = 1;
            else
                *id = 2;
            break;
        }
    }
    free(cand);
    return res;
}



void cont_fill_gasf(int *e, int f, int Q, fajor_sion3D * QUAD)
{
    int *cand, i, j, nd[4], n_a, n_b, nx, n1, n2, q;
    cand = (int *) malloc(7 * sizeof(int));
    for (i = 0; i < 6; i++)
        cand[i] = e[i];
    cand[6] = f;
    nd[0] = QUAD->elem[Q].frvrt;
    nd[1] = QUAD->elem[Q].scvrt;
    nd[2] = QUAD->elem[Q].thvrt;
    nd[3] = QUAD->elem[Q].ftvrt;
    for (i = 0; i < 4; i++) {
        n_a = nd[i];
        nx = i + 1;
        if (nx == 4)
            nx = 0;
        n_b = nd[nx];
        q = -1;
        for (j = 0; j < 7; j++) {
            n1 = QUAD->kt[cand[j]].frvrt;
            n2 = QUAD->kt[cand[j]].scvrt;
            if (((n1 == n_a) && (n2 == n_b)) || ((n1 == n_b) && (n2 == n_a))) {
                q = cand[j];
                break;
            }
        }
        if (q == -1) {
            fprintf(tmpout, "Flip: unable to find local edge\n");
            exit(0);
        }
        if (i == 0)
            QUAD->elem[Q].frkt = q;
        if (i == 1)
            QUAD->elem[Q].sckt = q;
        if (i == 2)
            QUAD->elem[Q].trkt = q;
        if (i == 3)
            QUAD->elem[Q].ftkt = q;
    }
    free(cand);
}



void jigk_flip_qirh(fajor_sion3D * quad, int *e, int f, int *nd, int p, int q)
{
    int E, F, i, z_e[4], z_f[4], nx, n_a, n_b, id;
    int cur_inc;
    fprintf(tmpout, "FLIPPING f=%d\n", f);
    E = quad->kt[f].frent;
    F = quad->kt[f].scent;

    z_e[0] = p;
    for (i = 1; i < 4; i++) {
        nx = nigw_peri_neql(z_e[i - 1]);
        z_e[i] = nx;
    }
    z_f[0] = q;
    for (i = 1; i < 4; i++) {
        nx = nigw_peri_neql(z_f[i - 1]);
        z_f[i] = nx;
    }
    for (i = 0; i < 4; i++) {
        bugq_assi_sokm(nd[z_e[i]], E, i + 1, quad);
        bugq_assi_sokm(nd[z_f[i]], F, i + 1, quad);
    }
    quad->kt[f].frvrt = nd[p];
    quad->kt[f].scvrt = nd[q];

    cont_fill_gasf(e, f, E, quad);
    cont_fill_gasf(e, f, F, quad);

    for (i = 0; i < 6; i++) {
        n_a = quad->kt[e[i]].frvrt;
        n_b = quad->kt[e[i]].scvrt;
        vogj_edge_tocf(n_a, n_b, *quad, E, F, &id);
        if (id == 1) {
            cur_inc = quad->kt[e[i]].frent;
            if (cur_inc == F)
                quad->kt[e[i]].frent = E;
            cur_inc = quad->kt[e[i]].scent;
            if (cur_inc == F)
                quad->kt[e[i]].scent = E;
        }
        if (id == 2) {
            cur_inc = quad->kt[e[i]].frent;
            if (cur_inc == E)
                quad->kt[e[i]].frent = F;
            cur_inc = quad->kt[e[i]].scent;
            if (cur_inc == E)
                quad->kt[e[i]].scent = F;
        }
    }
    gows_veri_jiqw(*quad);
}



double dugw_pros_jizk(float_curve * FC, fajor_sion3D quad, int *corresp, manif_tl msh, vect3D * nrm, int f, int *nd, int *e, int p, int q)
{
    int i, j, nx, w[4], E, F, dummy, h, cr, d;
    double q1, q2, res;
    point *cor, *mid, M;
    vect3D *N;

    nopr_best_jiwt(FC[f], &M);
    cor = (point *) malloc(4 * sizeof(point));
    mid = (point *) malloc(4 * sizeof(point));
    N = (vect3D *) malloc(4 * sizeof(vect3D));
    E = quad.kt[f].frent;
    F = quad.kt[f].scent;

    for (j = 0; j < 2; j++) {
        if (j == 0)
            cr = p;
        if (j == 1)
            cr = q;
        w[0] = nd[cr];
        for (i = 1; i < 4; i++) {
            nx = nigw_peri_neql(cr);
            cr = nx;
            w[i] = nd[cr];
        }
        for (i = 0; i < 4; i++) {
            getf_find_rogc_todj(quad.knot[w[i]], &cor[i]);
            d = corresp[w[i]];
            getf_find_rogc_todj(nrm[d], &N[i]);
        }
        for (i = 0; i < 4; i++) {
            nx = i + 1;
            if (nx == 4)
                nx = 0;
            h = vogj_edge_tocf(w[i], w[nx], quad, E, F, &dummy);
            if (h != -1)
                nopr_best_jiwt(FC[h], &mid[i]);
            else
                getf_find_rogc_todj(M, &mid[i]);
        }
        if (j == 0)
            q1 = lurq_qual_wotl(cor, mid, N, 1);
        if (j == 1)
            q2 = lurq_qual_wotl(cor, mid, N, 1);
    }
    free(cor);
    free(mid);
    free(N);
    if (q1 < q2)
        res = q1;
    else
        res = q2;
    return res;
}


void gicj_surr_qerd(fajor_sion3D quad, int f, int *nd)
{
    int i, nx, E, F, temp;
    E = quad.kt[f].frent;
    temp = cuzh_next_kadm(quad, E, quad.kt[f].frvrt);
    if (temp == quad.kt[f].scvrt)
        nd[0] = quad.kt[f].scvrt;
    else
        nd[0] = quad.kt[f].frvrt;
    for (i = 1; i < 4; i++) {
        nx = cuzh_next_kadm(quad, E, nd[i - 1]);
        nd[i] = nx;
    }
    F = quad.kt[f].scent;
    for (i = 4; i < 6; i++) {
        nx = cuzh_next_kadm(quad, F, nd[i - 1]);
        nd[i] = nx;
    }
}


void dekl_find_bofd(fajor_sion3D * quad, int *corresp, manif_tl msh, vect3D * nrm, float_curve * FC, int f, int *e, int *nd, double *LRG, double *CUR, int *P_SW, int *Q_SW)
{
    int i, j, nx, p, q, p_sw, q_sw, E, F, dummy;
    double lrg, cur, ql;
    gicj_surr_qerd(*quad, f, nd);
    E = quad->kt[f].frent;
    F = quad->kt[f].scent;
    for (i = 0; i < 6; i++) {
        nx = i + 1;
        if (nx == 6)
            nx = 0;
        e[i] = vogj_edge_tocf(nd[i], nd[nx], *quad, E, F, &dummy);
        if (e[i] == -1) {
            fprintf(tmpout, "separator: f[%d]=[%d,%d]\n", f, quad->kt[f].frvrt, quad->kt[f].scvrt);
            jans_disp_nudj(*quad, E);
            jans_disp_nudj(*quad, F);
            for (j = 0; j < 6; j++)
                fprintf(tmpout, "nd[%d]=%d\n", j, nd[j]);
            fprintf(tmpout, "Unable to find hex bound kt  [%d,%d]\n", nd[i], nd[nx]);
            exit(0);
        }
    }
    cur = dugw_pros_jizk(FC, *quad, corresp, msh, nrm, f, nd, e, 0, 3);
    lrg = -LARGE_NUMBER;
    for (p = 1; p <= 2; p++) {
        q = p + 3;
        ql = dugw_pros_jizk(FC, *quad, corresp, msh, nrm, f, nd, e, p, q);
        if (ql > lrg) {
            lrg = ql;
            p_sw = p;
            q_sw = q;
        }
    }
    *LRG = lrg;
    *CUR = cur;
    *P_SW = p_sw;
    *Q_SW = q_sw;
}


int cfh = 0;


void seph_cond_guzw(float_curve * FC, fajor_sion3D * quad, int f, prat_main GR, manif_tl msh, vect3D * nrm, int *corresp, int max_nb_smooth)
{
    int *e, *nd, p_sw, q_sw, nb_inter = 5;
    int N1, N2, n1, n2, level = 0, nb_stat, *type, *underl;
    double cur, lrg;
    point *omega;
    e = (int *) malloc(6 * sizeof(int));
    nd = (int *) malloc(6 * sizeof(int));
    dekl_find_bofd(quad, corresp, msh, nrm, FC, f, e, nd, &lrg, &cur, &p_sw, &q_sw);
    if (lrg > cur) {
        jigk_flip_qirh(quad, e, f, nd, p_sw, q_sw);
        N1 = quad->kt[f].frvrt;
        N2 = quad->kt[f].scvrt;
        n1 = corresp[N1];
        n2 = corresp[N2];
        omega = (point *) malloc(max_nb_smooth * sizeof(point));
        type = (int *) malloc(max_nb_smooth * sizeof(int));
        underl = (int *) malloc(max_nb_smooth * sizeof(int));
        nb_stat = folc_gene_kost(n1, n2, level, nb_inter, GR, msh, omega, type, underl, max_nb_smooth);
        with_conv_davf(msh, GR, omega, type, underl, nb_stat, &FC[f]);
        if (FC[f].st_grs >= max_nb_smooth) {
            fprintf(tmpout, "2-Max number of stations [%d]\n", max_nb_smooth);
            exit(0);
        }
        zekl_redu_govr(msh, &FC[f]);
        free(omega);
        free(type);
        free(underl);

        cfh++;
    }
    free(nd);
    free(e);
}


void nuws_impr_gudk(float_curve * FC, fajor_sion3D * quad, prat_main GR, manif_tl msh, vect3D * nrm, int *corresp, int max_nb_smooth)
{
    int f, ned, E, F, ts, dummy;
    fprintf(tmpout, "IMPROVE FLIP\n");
    gows_veri_jiqw(*quad);
    ned = quad->k_grs;
    for (f = 0; f < ned; f++) {
        E = quad->kt[f].frent;
        F = quad->kt[f].scent;
        ts = cojs_shar_sejm(*quad, E, F, &dummy);
        if ((ts == 0) && (E != F))
            seph_cond_guzw(FC, quad, f, GR, msh, nrm, corresp, max_nb_smooth);
    }
}



void lutq_find_pobh_mugq(int w, fajor_sion3D * quad, float_curve * FC, prat_main GR, manif_tl msh, int *corresp, vect3D * nrm, int max_nb_smooth)
{
    int f[4], i, j, q = -1, *e, *nd, p_sw, q_sw, P_SW, Q_SW, *E, *ND, e1, e2, ts;
    int N1, N2, n1, n2, *type, *underl, nb_stat, level = 0, nb_inter = 5, obs;
    double cur, lrg, LRG;
    point *omega;
    f[0] = quad->elem[w].frkt;
    f[1] = quad->elem[w].sckt;
    f[2] = quad->elem[w].trkt;
    f[3] = quad->elem[w].ftkt;
    e = (int *) malloc(6 * sizeof(int));
    nd = (int *) malloc(6 * sizeof(int));
    E = (int *) malloc(6 * sizeof(int));
    ND = (int *) malloc(6 * sizeof(int));
    LRG = -LARGE_NUMBER;
    for (i = 0; i < 4; i++) {
        e1 = quad->kt[f[i]].frent;
        e2 = quad->kt[f[i]].scent;
        ts = cojs_shar_sejm(*quad, e1, e2, &obs);
        if (ts == 0) {
            dekl_find_bofd(quad, corresp, msh, nrm, FC, f[i], e, nd, &lrg, &cur, &p_sw, &q_sw);
            if ((lrg > cur) && (lrg > LRG)) {
                q = i;
                LRG = lrg;
                P_SW = p_sw;
                Q_SW = q_sw;
                for (j = 0; j < 6; j++) {
                    E[j] = e[j];
                    ND[j] = nd[j];
                }
            }
        }
    }
    if (q != -1) {
        jigk_flip_qirh(quad, E, f[q], ND, P_SW, Q_SW);
        N1 = quad->kt[f[q]].frvrt;
        N2 = quad->kt[f[q]].scvrt;
        n1 = corresp[N1];
        n2 = corresp[N2];
        omega = (point *) malloc(max_nb_smooth * sizeof(point));
        type = (int *) malloc(max_nb_smooth * sizeof(int));
        underl = (int *) malloc(max_nb_smooth * sizeof(int));
        nb_stat = folc_gene_kost(n1, n2, level, nb_inter, GR, msh, omega, type, underl, max_nb_smooth);
        with_conv_davf(msh, GR, omega, type, underl, nb_stat, &FC[f[q]]);
        if (FC[f[q]].st_grs >= max_nb_smooth) {
            fprintf(tmpout, "3-Max number of stations [%d]\n", max_nb_smooth);
            exit(0);
        }
        zekl_redu_govr(msh, &FC[f[q]]);
        free(omega);
        free(type);
        free(underl);
    }
    free(nd);
    free(e);
    free(ND);
    free(E);
}


void vizr_find_jumr_nazp(prat_main GR, fajor_sion3D * quad, float_curve * FC, manif_tl msh, int *corresp, vect3D * nrm, int max_nb_smooth)
{
    int *E, i, k;
    double Q;
    E = (int *) malloc(quad->e_grs * sizeof(int));
    k = 0;
    for (i = 0; i < quad->e_grs; i++) {
        Q = lijf_loca_necf(FC, *quad, corresp, nrm, i);
        if (Q < 0.0) {
            fprintf(tmpout, "Element=%d  quality=%f\n", i, Q);
            E[k] = i;
            lutq_find_pobh_mugq(E[k], quad, FC, GR, msh, corresp, nrm, max_nb_smooth);
            k++;
        }
    }
    free(E);
}
