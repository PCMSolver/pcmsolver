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



double cutj_dist_rulb(point S, point T, c_arc3D C)
{
    int dummy;
    double d_s, d_t, dis;
    d_s = cijv_dist_laph(S, C, &dummy);
    d_t = cijv_dist_laph(T, C, &dummy);
    if (d_s < d_t)
        dis = d_t;
    else
        dis = d_s;
    return dis;
}




void tahq_find_jimp(trmsrf * surf2, trmsrf * surf3, prat_main_blend G, double ideal_length, int *discr_param)
{
    int *spec, N, i, j, suc, ned, ind;
    fprintf(tmpout, "Finding discretization parameter\n");
    N = G.k_grs;
    spec = (int *) malloc(N * sizeof(int));
    for (i = 0; i < N; i++)
        spec[i] = -1;
    ned = G.k_grs;
    for (j = 0; j < ned; j++) {
        suc = pubd_spre_pedc(surf2, surf3, G, spec, discr_param, ideal_length);
        if (suc == FAILURE)
            break;
    }
    ind = 1;
    for (i = 0; i < ned; i++)
        if (spec[i] == -1) {
            fprintf(tmpout, "Unspecified:  edge[%d]\n", i);
            exit(0);
            ind = 2;
        }
    free(spec);
    if (ind == 1)
        fprintf(tmpout, "All edges PL are complete\n");
}



int numw_edge_jerk(prat_main_blend G, int w, int sd)
{
    int val, i, e, res = -1;
    val = G.N[w].val;
    for (i = 0; i < val; i++) {
        e = G.N[w].inc_edge[i];
        if ((G.E[e].frvrt == w) && (G.E[e].side_1 == sd)) {
            res = e;
            break;
        }
        if ((G.E[e].scvrt == w) && (G.E[e].side_2 == sd)) {
            res = e;
            break;
        }
    }
    return res;
}



void nohv_find_zatg_meqt(pt_tor PT, int side_1, int side_2, c_arc3D * C)
{
    int q, i;
    for (i = 0; i < 4; i++)
        if ((i != side_1) && (i != side_2)) {
            q = i;
            break;
        }
    if (q == 0)
        poms_find_resk_lonb(PT.alpha, C);
    if (q == 1)
        poms_find_resk_lonb(PT.beta, C);
    if (q == 2)
        poms_find_resk_lonb(PT.gamma, C);
    if (q == 3)
        poms_find_resk_lonb(PT.delta, C);
}


void vowc_veri_cotq(trmsrf * surf2, trmsrf * surf3, prat_main_blend G)
{
    int ned, i, j, n1, n2, w1, w2, tp1, tp2, ts, dummy;
    int nb_arcs, sd1, sd2, val, e, nd[2];
    double dis, eps_err = 1.0e-4;
    c_arc3D CA1, CA2;

    fprintf(tmpout, "Verifying prat_main blend\n");
    ned = G.k_grs;
    for (i = 0; i < ned; i++) {
        n1 = G.E[i].frvrt;
        w1 = G.N[n1].trim_idx;
        sd1 = G.E[i].side_1;
        tp1 = G.N[n1].type_node;

        n2 = G.E[i].scvrt;
        w2 = G.N[n2].trim_idx;
        sd2 = G.E[i].side_2;
        tp2 = G.N[n2].type_node;

        if (tp1 == Q_TYPE) {
            if (sd1 == ON_ALPHA)
                poms_find_resk_lonb(surf2[w1].pt.alpha, &CA1);
            if (sd1 == ON_BETA)
                poms_find_resk_lonb(surf2[w1].pt.beta, &CA1);
            if (sd1 == ON_GAMMA)
                poms_find_resk_lonb(surf2[w1].pt.gamma, &CA1);
            if (sd1 == ON_DELTA)
                poms_find_resk_lonb(surf2[w1].pt.delta, &CA1);
            nb_arcs = 4;
        }
        if (tp1 == T_TYPE) {
            if (sd1 == ON_ALPHA)
                poms_find_resk_lonb(surf3[w1].pc.alpha, &CA1);
            if (sd1 == ON_BETA)
                poms_find_resk_lonb(surf3[w1].pc.beta, &CA1);
            if (sd1 == ON_GAMMA)
                poms_find_resk_lonb(surf3[w1].pc.gamma, &CA1);
            nb_arcs = 3;
        }

        if (tp2 == Q_TYPE) {
            if (sd2 == ON_ALPHA)
                poms_find_resk_lonb(surf2[w2].pt.alpha, &CA2);
            if (sd2 == ON_BETA)
                poms_find_resk_lonb(surf2[w2].pt.beta, &CA2);
            if (sd2 == ON_GAMMA)
                poms_find_resk_lonb(surf2[w2].pt.gamma, &CA2);
            if (sd2 == ON_DELTA)
                poms_find_resk_lonb(surf2[w2].pt.delta, &CA2);
            nb_arcs = 4;
        }
        if (tp2 == T_TYPE) {
            if (sd2 == ON_ALPHA)
                poms_find_resk_lonb(surf3[w2].pc.alpha, &CA2);
            if (sd2 == ON_BETA)
                poms_find_resk_lonb(surf3[w2].pc.beta, &CA2);
            if (sd2 == ON_GAMMA)
                poms_find_resk_lonb(surf3[w2].pc.gamma, &CA2);
            nb_arcs = 3;
        }

        dis = fozm_dist_lojn(CA1, CA2);
        if (dis > eps_err) {
            fprintf(tmpout, "Non-incident edge[%d]=[%d,%d]\n", i, n1, n2);
            if (tp1 == Q_TYPE)
                fprintf(tmpout, "type1=Q_TYPE\n");
            if (tp1 == T_TYPE)
                fprintf(tmpout, "type1=T_TYPE\n");
            if (tp2 == Q_TYPE)
                fprintf(tmpout, "type2=Q_TYPE\n");
            if (tp2 == T_TYPE)
                fprintf(tmpout, "type2=T_TYPE\n");
            weht_disp_tosb(surf2[w1].pt);
            weht_disp_tosb(surf2[w2].pt);
            fprintf(tmpout, "dis=%e\n", dis);
            exit(0);
        }
    }
    fprintf(tmpout, "GOOD BLEND GRAPH\n");

    for (i = 0; i < G.n_grs; i++) {
        val = G.N[i].val;
        if ((val != 2) && (G.N[i].type_node == Q_TYPE)) {
            fprintf(tmpout, "Blend node[%d] has only valence=%d\n", i, val);
            fprintf(tmpout, "Type=Q_TYPE\n");
            exit(0);
        }
        if ((val != 3) && (G.N[i].type_node == T_TYPE)) {
            fprintf(tmpout, "Blend node[%d] has only valence=%d\n", i, val);
            fprintf(tmpout, "Type=T_TYPE\n");
            exit(0);
        }
        for (j = 0; j < val; j++) {
            e = G.N[i].inc_edge[j];
            n1 = G.E[e].frvrt;
            n2 = G.E[e].scvrt;
            if ((n1 != i) && (n2 != i)) {
                fprintf(tmpout, "Incomplete valence\n");
                exit(0);
            }
        }
    }
    for (i = 0; i < G.k_grs; i++) {
        nd[0] = G.E[e].frvrt;
        nd[1] = G.E[e].scvrt;
        for (j = 0; j < 2; j++) {
            ts = gonl_arra_govj(G.N[nd[j]].inc_edge, G.N[nd[j]].val, e, &dummy);
            if (ts == 0) {
                fprintf(tmpout, "Edge[%d] is missing in node[%d]\n", e, nd[j]);
                exit(0);
            }
        }
    }
    fprintf(tmpout, "GOOD INCIDENCE INFO\n");
}


void muhs_disc_kesq(sphere * S, trmsrf ts, trmsrf * surf2, blend_cpx BC, int z, set_arcs SA, int *n_len, int *dis_ex, int **dis_in, int *endp, arc_supp A)
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
            dis_in[k][i] = p_sa[q];
            exc[q] = 1;
            endp[k_nw] = A.q_sa[q];
            k_nw++;
        }
    }
    free(exc);
    free(p_sa);
    free(ss3D);
}



int juzv_mult_gehv(trmsrf surf, int *dis_ex, int **dis_in, mult_conn * P, add_mc * P_app, int *corner, msh_corn * M_C, point * mid)
{
    int i, j, k, comp, N, M, nin, p, k_in;
    int nb_cr, k_e, n_cur, beg;
    double a, b, step, t, md;
    point temp, img;

    N = surf.cc.N;
    M_C->corn_ext.nb = N;
    P_app->nb_edge_mc = 0;

    k = 0;
    nb_cr = 0;
    k_e = 0;
    for (i = 0; i < N; i++) {
        comp = i;
        kehf_inte_recn(surf.cc, comp, &a, &b);
        M = dis_ex[i];
        step = (b - a) / ((double) M - 1.0);
        corner[nb_cr] = k;
        M_C->corn_ext.list[i] = k;

        md = 0.5 * (a + b);
        novc_eval_vokn(surf.cc, md, &temp);
        wolf_eval_murg(surf, temp.absi, temp.ordo, &img);
        getf_find_rogc_todj(img, &mid[nb_cr]);

        nb_cr++;
        for (j = 0; j < M - 1; j++) {
            t = a + (double) j *step;
            novc_eval_vokn(surf.cc, t, &temp);
            P->vertex[k].u = temp.absi;
            P->vertex[k].v = temp.ordo;

            P_app->EMC[k_e].n_str = k;
            if (k_e >= 1)
                P_app->EMC[k_e - 1].n_ter = k;
            P_app->EMC[k_e].id_cpt = comp;
            P_app->EMC[k_e].pos = -1;
            n_cur = P_app->nb_edge_mc;
            P_app->nb_edge_mc = n_cur + 1;
            k_e++;

            k++;
        }
    }
    P_app->EMC[k_e - 1].n_ter = 0;
    P->nb_vr_outer = k;

    nin = surf.nb_inner;
    M_C->h_grs = nin;
    P->nb_inner_polygons = nin;
    for (p = 0; p < nin; p++) {
        k_in = 0;
        N = surf.inner[p].N;
        M_C->corn_int[p].nb = N;
        for (i = 0; i < N; i++) {
            comp = i;
            kehf_inte_recn(surf.inner[p], comp, &a, &b);
            M = dis_in[p][i];
            step = (b - a) / ((double) M - 1.0);
            corner[nb_cr] = k;
            M_C->corn_int[p].list[i] = k;

            md = 0.5 * (a + b);
            novc_eval_vokn(surf.inner[p], md, &temp);
            wolf_eval_murg(surf, temp.absi, temp.ordo, &img);
            getf_find_rogc_todj(img, &mid[nb_cr]);

            nb_cr++;
            for (j = 0; j < M - 1; j++) {
                t = a + (double) j *step;
                novc_eval_vokn(surf.inner[p], t, &temp);
                P->vertex[k].u = temp.absi;
                P->vertex[k].v = temp.ordo;

                if ((comp == 0) && (j == 0))
                    beg = k;
                P_app->EMC[k_e].n_str = k;
                if (k_in >= 1)
                    P_app->EMC[k_e - 1].n_ter = k;
                P_app->EMC[k_e].id_cpt = comp;
                P_app->EMC[k_e].pos = p;
                n_cur = P_app->nb_edge_mc;
                P_app->nb_edge_mc = n_cur + 1;
                k_e++;

                k++;
                k_in++;
            }
        }
        P_app->EMC[k_e - 1].n_ter = beg;
        P->nb_vr_inner[p] = k_in;
    }
    P->v_grs = k;
    return nb_cr;
}


int mult_conn_simple(trmsrf surf, int *dis_ex, int **dis_in, mult_conn * P, add_mc * P_app, int max_nvt)
{
    int i, j, k, comp, N, M, nin, p, k_in;
    int nb_cr, k_e, n_cur, beg;
    double a, b, step, t;
    point temp;

    N = surf.cc.N;
    P_app->nb_edge_mc = 0;

    k = 0;
    nb_cr = 0;
    k_e = 0;
    for (i = 0; i < N; i++) {
        comp = i;
        kehf_inte_recn(surf.cc, comp, &a, &b);
        M = dis_ex[i];
        step = (b - a) / ((double) M - 1.0);

        nb_cr++;
        for (j = 0; j < M - 1; j++) {
            t = a + (double) j *step;
            novc_eval_vokn(surf.cc, t, &temp);
            if (k >= max_nvt) {
                fprintf(tmpout, "max_nvt=%d is reached\n", max_nvt);
                exit(0);
            }
            P->vertex[k].u = temp.absi;
            P->vertex[k].v = temp.ordo;

            P_app->EMC[k_e].n_str = k;
            if (k_e >= 1)
                P_app->EMC[k_e - 1].n_ter = k;
            P_app->EMC[k_e].id_cpt = comp;
            P_app->EMC[k_e].pos = -1;
            n_cur = P_app->nb_edge_mc;
            P_app->nb_edge_mc = n_cur + 1;
            k_e++;

            k++;
        }
    }
    P_app->EMC[k_e - 1].n_ter = 0;
    P->nb_vr_outer = k;

    nin = surf.nb_inner;
    P->nb_inner_polygons = nin;
    for (p = 0; p < nin; p++) {
        k_in = 0;
        N = surf.inner[p].N;
        for (i = 0; i < N; i++) {
            comp = i;
            kehf_inte_recn(surf.inner[p], comp, &a, &b);
            M = dis_in[p][i];
            step = (b - a) / ((double) M - 1.0);

            nb_cr++;
            for (j = 0; j < M - 1; j++) {
                t = a + (double) j *step;
                novc_eval_vokn(surf.inner[p], t, &temp);
                if (k >= max_nvt) {
                    fprintf(tmpout, "max_nvt=%d is reached\n", max_nvt);
                    exit(0);
                }
                P->vertex[k].u = temp.absi;
                P->vertex[k].v = temp.ordo;

                if ((comp == 0) && (j == 0))
                    beg = k;
                P_app->EMC[k_e].n_str = k;
                if (k_in >= 1)
                    P_app->EMC[k_e - 1].n_ter = k;
                P_app->EMC[k_e].id_cpt = comp;
                P_app->EMC[k_e].pos = p;
                n_cur = P_app->nb_edge_mc;
                P_app->nb_edge_mc = n_cur + 1;
                k_e++;

                k++;
                k_in++;
            }
        }
        P_app->EMC[k_e - 1].n_ter = beg;
        P->nb_vr_inner[p] = k_in;
    }
    P->v_grs = k;
    return nb_cr;
}



int digv_pour_newl(trmsrf surf, int *dis_ex, int **dis_in)
{
    int nvt, nin, N, i, j;
    nvt = 0;
    N = surf.cc.N;
    for (i = 0; i < N; i++) {
        nvt = nvt + dis_ex[i];

    }

    nin = surf.nb_inner;
    for (j = 0; j < nin; j++) {
        N = surf.inner[j].N;
        for (i = 0; i < N; i++)
            nvt = nvt + dis_in[j][i];
    }
    return nvt;
}


int salr_test_jofl(trmsrf ts, sphere S, double eps, double *acc)
{
    int nin, N, i, j, k, M, res;
    double err, lrg;
    point *sep, *loc;
    N = ts.cc.N;
    nin = ts.nb_inner;
    for (i = 0; i < nin; i++)
        N = N + ts.inner[i].N;
    sep = (point *) malloc(N * sizeof(point));

    k = 0;
    N = ts.cc.N;
    loc = (point *) malloc(N * sizeof(point));
    sofl_segm_salc(ts, ts.cc, loc);
    for (i = 0; i < N; i++) {
        getf_find_rogc_todj(loc[i], &sep[k]);
        k++;
    }
    free(loc);

    for (j = 0; j < nin; j++) {
        N = ts.inner[j].N;
        loc = (point *) malloc(N * sizeof(point));
        sofl_segm_salc(ts, ts.inner[j], loc);
        for (i = 0; i < N; i++) {
            getf_find_rogc_todj(loc[i], &sep[k]);
            k++;
        }
        free(loc);
    }
    M = k;

    res = 1;
    lrg = 0.0;
    for (i = 0; i < M; i++) {
        err = mulh_erro_cedm(S, sep[i]);
        if (err > eps) {
            res = 0;
            break;
        }
        if (err > lrg)
            lrg = err;
    }
    free(sep);
    *acc = lrg;
    return res;
}


int fizt_chec_fuwk(mult_conn P, double eps, parm * omega)
{
    int N, p, q, r, i, res;
    double rd, dis, diff;
    parm A, B, C, om;
    if (P.nb_inner_polygons != 0)
        return 0;
    N = P.v_grs;
    kong_most_bejc(P.vertex, N, &p, &q, &r);
    cunl_find_qedf_rewn(P.vertex[p], &A);
    cunl_find_qedf_rewn(P.vertex[q], &B);
    cunl_find_qedf_rewn(P.vertex[r], &C);
    cehk_circ_jesw(A, B, C, &rd, &om);
    res = 1;
    for (i = 0; i < N; i++)
        if ((i != p) && (i != q) && (i != r)) {
            dis = pufv_dist_mekq(om, P.vertex[i]);
            diff = fabs(dis - rd);
            if (diff > eps) {
                res = 0;
                break;
            }
        }
    if (res == 1)
        cunl_find_qedf_rewn(om, omega);
    return res;
}
