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
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "sas.h"
#include "smooth.h"


int lucw_find_feph_qofd(manif_tl msh, prat_main GR, int s, int *ngb)
{
    int i, e, nb, val, ts, dummy, nd1, nd2, m;
    nb = 0;
    val = GR.dgr[s];
    for (i = 0; i < val; i++) {
        e = GR.incd[s][i];
        nd1 = msh.kt[e].frvrt;
        nd2 = msh.kt[e].scvrt;
        if (nd1 == s)
            m = nd2;
        if (nd2 == s)
            m = nd1;
        ts = gonl_arra_govj(ngb, nb, m, &dummy);
        if (ts == 0) {
            ngb[nb] = m;
            nb++;
        }
    }
    return nb;
}



int zacr_find_netw_wusf(int s, prat_main GR, manif_tl msh, int level, int *ngb, int max_size)
{
    int nb, i, j, k, *temp, val, beg, z, nb_loc, ts, dummy, nb_old;
    nb = lucw_find_feph_qofd(msh, GR, s, ngb);
    if (nb >= max_size) {
        fprintf(tmpout, "max_size is obtained\n");
        exit(0);
    }
    beg = 0;
    if (level >= 1)
        for (i = 1; i <= level; i++) {
            nb_old = nb;
            for (j = beg; j < nb_old; j++) {
                z = ngb[j];
                val = GR.dgr[z];
                temp = (int *) malloc(val * sizeof(int));
                nb_loc = lucw_find_feph_qofd(msh, GR, z, temp);
                for (k = 0; k < nb_loc; k++) {
                    ts = gonl_arra_govj(ngb, nb, temp[k], &dummy);
                    if (ts == 0) {
                        if (nb >= max_size) {
                            fprintf(tmpout, "max_size is obtained\n");
                            exit(0);
                        }
                        ngb[nb] = temp[k];
                        nb++;
                    }
                }
                free(temp);
            }
            beg = nb_old;
            nb_old = nb;
        }
    return nb;
}



int cahq_idea_fewg(fajor_sion3D quad, int *corresp, prat_main GR, manif_tl msh, int lev, int nd_a, int nd_b)
{
    int max_size = 100, *cand, nb_cand, max_cand = 1000, i, z;
    int id, m_a, m_b, *ngb, N, ts, dummy;
    double sml, dis;
    point X;

    X.absi = 0.5 * (quad.knot[nd_a].absi + quad.knot[nd_b].absi);
    X.ordo = 0.5 * (quad.knot[nd_a].ordo + quad.knot[nd_b].ordo);
    X.cote = 0.5 * (quad.knot[nd_a].cote + quad.knot[nd_b].cote);

    m_a = corresp[nd_a];
    m_b = corresp[nd_b];
    cand = (int *) malloc(max_cand * sizeof(int));
    ngb = (int *) malloc(max_size * sizeof(int));
    N = zacr_find_netw_wusf(m_a, GR, msh, lev, ngb, max_size);
    for (i = 0; i < N; i++)
        cand[i] = ngb[i];
    nb_cand = N;
    N = zacr_find_netw_wusf(m_b, GR, msh, lev, ngb, max_size);
    for (i = 0; i < N; i++) {
        ts = gonl_arra_govj(cand, nb_cand, ngb[i], &dummy);
        if (ts == 0) {
            cand[nb_cand] = ngb[i];
            nb_cand++;
        }
    }
    free(ngb);

    sml = LARGE_NUMBER;
    for (i = 0; i < nb_cand; i++) {
        z = cand[i];
        dis = wodt_dist_gilq(X, msh.knot[z]);
        if (dis < sml) {
            sml = dis;
            id = z;
        }
    }
    free(cand);
    return id;
}


void ripn_find_reln(fajor_sion3D quad, int F, int nd_a, int nd_b, int *m_1, int *m_2, int *e_a, int *e_b, int *f_a, int *f_b)
{
    int M1 = -1, M2 = -1, nd[4], e[4], i, gr;
    int E_a, E_b, F_a, F_b, n1, n2;
    nd[0] = quad.elem[F].frvrt;
    nd[1] = quad.elem[F].scvrt;
    nd[2] = quad.elem[F].thvrt;
    nd[3] = quad.elem[F].ftvrt;
    for (i = 0; i < 4; i++)
        if ((nd[i] != nd_a) && (nd[i] != nd_b)) {
            gr = M1;
            if (gr == -1)
                M1 = nd[i];
            else
                M2 = nd[i];
            if ((M1 != -1) && (M2 != -1))
                break;
        }
    *m_1 = M1;
    *m_2 = M2;

    e[0] = quad.elem[F].frkt;
    e[1] = quad.elem[F].sckt;
    e[2] = quad.elem[F].trkt;
    e[3] = quad.elem[F].ftkt;
    for (i = 0; i < 4; i++) {
        n1 = quad.kt[e[i]].frvrt;
        n2 = quad.kt[e[i]].scvrt;
        if (((n1 == M1) && (n2 == nd_a)) || ((n2 == M1) && (n1 == nd_a)))
            E_a = e[i];
        if (((n1 == M1) && (n2 == nd_b)) || ((n2 == M1) && (n1 == nd_b)))
            F_a = e[i];
        if (((n1 == M2) && (n2 == nd_a)) || ((n2 == M2) && (n1 == nd_a)))
            E_b = e[i];
        if (((n1 == M2) && (n2 == nd_b)) || ((n2 == M2) && (n1 == nd_b)))
            F_b = e[i];
    }
    *e_a = E_a;
    *e_b = E_b;
    *f_a = F_a;
    *f_b = F_b;
}


int minf_invo_jawc(fajor_sion3D quad, int nd_a, int nd_b, int e_a, int f_a, int e_b, int f_b, int *inv_ed, int max_inv)
{
    int ned, i, n1, n2, N;
    ned = quad.k_grs;
    N = 0;
    for (i = 0; i < ned; i++)
        if ((i != e_a) && (i != e_b) && (i != f_a) && (i != f_b)) {
            n1 = quad.kt[i].frvrt;
            n2 = quad.kt[i].scvrt;
            if ((n1 == nd_a) || (n2 == nd_a) || (n1 == nd_b) || (n2 == nd_b)) {
                if (N >= max_inv) {
                    fprintf(tmpout, "max_inv is reached\n");
                    exit(0);
                }
                inv_ed[N] = i;
                N++;
            }
        }
    return N;
}


void jeqg_skit_leqb(fajor_sion3D * quad, manif_tl msh, prat_main GR, int F, int nd_a, int nd_b, int id, int *corresp, float_curve * fc, int nb_inter, int max_nb_smooth)
{
    int m_1, m_2, e_a, e_b, f_a, f_b, i, k, inc_e, new_e_a, new_e_b;
    int n1, n2, n3, n4, e1, e2, e3, e4, q1, q2, ned, z, nb_stat;
    int obs, inc_f, *map_ed, *map_qd, *map_nd, level = 0, *type, *underl;
    int im_n1, im_n2, *inv_ed, max_inv = 50, n_inv, temp;
    kt_t *temp_ed;
    float_curve *temp_fc;
    efajor *temp_qd;
    point *temp_nd, *omega;

    ripn_find_reln(*quad, F, nd_a, nd_b, &m_1, &m_2, &e_a, &e_b, &f_a, &f_b);
    inv_ed = (int *) malloc((max_inv + 2) * sizeof(int));
    n_inv = minf_invo_jawc(*quad, nd_a, nd_b, e_a, f_a, e_b, f_b, inv_ed, max_inv);
    getf_find_rogc_todj(msh.knot[id], &quad->knot[nd_a]);

    for (i = 0; i < quad->e_grs; i++) {
        n1 = quad->elem[i].frvrt;
        if (n1 == nd_b)
            quad->elem[i].frvrt = nd_a;
        n2 = quad->elem[i].scvrt;
        if (n2 == nd_b)
            quad->elem[i].scvrt = nd_a;
        n3 = quad->elem[i].thvrt;
        if (n3 == nd_b)
            quad->elem[i].thvrt = nd_a;
        n4 = quad->elem[i].ftvrt;
        if (n4 == nd_b)
            quad->elem[i].ftvrt = nd_a;

        e1 = quad->elem[i].frkt;
        if (e1 == f_a)
            quad->elem[i].frkt = e_a;
        e2 = quad->elem[i].sckt;
        if (e2 == f_a)
            quad->elem[i].sckt = e_a;
        e3 = quad->elem[i].trkt;
        if (e3 == f_a)
            quad->elem[i].trkt = e_a;
        e4 = quad->elem[i].ftkt;
        if (e4 == f_a)
            quad->elem[i].ftkt = e_a;

        e1 = quad->elem[i].frkt;
        if (e1 == f_b)
            quad->elem[i].frkt = e_b;
        e2 = quad->elem[i].sckt;
        if (e2 == f_b)
            quad->elem[i].sckt = e_b;
        e3 = quad->elem[i].trkt;
        if (e3 == f_b)
            quad->elem[i].trkt = e_b;
        e4 = quad->elem[i].ftkt;
        if (e4 == f_b)
            quad->elem[i].ftkt = e_b;
    }

    if ((quad->kt[e_a].frent != F) && (quad->kt[e_a].scent != F)) {
        fprintf(tmpout, "F must be incident upon e_a\n");
        exit(0);
    }
    if ((quad->kt[f_a].frent != F) && (quad->kt[f_a].scent != F)) {
        fprintf(tmpout, "F must be incident upon f_a\n");
        exit(0);
    }
    if (quad->kt[e_a].frent == F)
        inc_e = quad->kt[e_a].scent;
    else
        inc_e = quad->kt[e_a].frent;
    if (quad->kt[f_a].frent == F)
        inc_f = quad->kt[f_a].scent;
    else
        inc_f = quad->kt[f_a].frent;
    quad->kt[e_a].frent = inc_e;
    quad->kt[e_a].scent = inc_f;

    if ((quad->kt[e_b].frent != F) && (quad->kt[e_b].scent != F)) {
        fprintf(tmpout, "F must be incident upon e_b\n");
        exit(0);
    }
    if ((quad->kt[f_b].frent != F) && (quad->kt[f_b].scent != F)) {
        fprintf(tmpout, "F must be incident upon f_b\n");
        exit(0);
    }
    if (quad->kt[e_b].frent == F)
        inc_e = quad->kt[e_b].scent;
    else
        inc_e = quad->kt[e_b].frent;
    if (quad->kt[f_b].frent == F)
        inc_f = quad->kt[f_b].scent;
    else
        inc_f = quad->kt[f_b].frent;
    quad->kt[e_b].frent = inc_e;
    quad->kt[e_b].scent = inc_f;

    for (i = 0; i < quad->k_grs; i++) {
        n1 = quad->kt[i].frvrt;
        if (n1 == nd_b)
            quad->kt[i].frvrt = nd_a;
        n2 = quad->kt[i].scvrt;
        if (n2 == nd_b)
            quad->kt[i].scvrt = nd_a;
    }

    q1 = f_a;
    q2 = f_b;
    ned = quad->k_grs;
    temp_fc = (float_curve *) malloc(ned * sizeof(float_curve));
    for (i = 0; i < ned; i++)
        vewk_allo_jovk(max_nb_smooth, &temp_fc[i]);
    temp_ed = (kt_t *) malloc(ned * sizeof(kt_t));
    map_ed = (int *) malloc(ned * sizeof(int));
    k = 0;
    for (i = 0; i < ned; i++)
        if ((i != q1) && (i != q2)) {
            if (fc[i].st_grs >= max_nb_smooth) {
                fprintf(tmpout, "max_nb_smooth is exceeded\n");
                exit(0);
            }
            rekf_find_kujn_rukg(fc[i], &temp_fc[k]);
            temp_ed[k].frvrt = quad->kt[i].frvrt;
            temp_ed[k].scvrt = quad->kt[i].scvrt;
            temp_ed[k].frent = quad->kt[i].frent;
            temp_ed[k].scent = quad->kt[i].scent;
            map_ed[i] = k;
            k++;
        }
    map_ed[q1] = -1;
    map_ed[q2] = -1;

    for (i = 0; i < k; i++) {
        rekf_find_kujn_rukg(temp_fc[i], &fc[i]);
        quad->kt[i].frvrt = temp_ed[i].frvrt;
        quad->kt[i].scvrt = temp_ed[i].scvrt;
        quad->kt[i].frent = temp_ed[i].frent;
        quad->kt[i].scent = temp_ed[i].scent;
    }
    quad->k_grs = k;
    for (i = 0; i < ned; i++)
        lohm_dest_nosr(&temp_fc[i]);
    free(temp_fc);
    free(temp_ed);
    for (i = 0; i < quad->e_grs; i++) {
        e1 = quad->elem[i].frkt;
        quad->elem[i].frkt = map_ed[e1];
        e2 = quad->elem[i].sckt;
        quad->elem[i].sckt = map_ed[e2];
        e3 = quad->elem[i].trkt;
        quad->elem[i].trkt = map_ed[e3];
        e4 = quad->elem[i].ftkt;
        quad->elem[i].ftkt = map_ed[e4];
    }
    new_e_a = map_ed[e_a];
    new_e_b = map_ed[e_b];
    for (i = 0; i < n_inv; i++) {
        temp = map_ed[inv_ed[i]];
        inv_ed[i] = temp;
    }
    free(map_ed);

    temp_qd = (efajor *) malloc(quad->e_grs * sizeof(efajor));
    map_qd = (int *) malloc(quad->e_grs * sizeof(int));
    k = 0;
    for (i = 0; i < quad->e_grs; i++)
        if (i != F) {
            fupj_find_numk_jobd(quad->elem[i], &temp_qd[k]);
            map_qd[i] = k;
            k++;
        }

    for (i = 0; i < k; i++)
        fupj_find_numk_jobd(temp_qd[i], &quad->elem[i]);
    quad->e_grs = k;
    free(temp_qd);
    for (i = 0; i < quad->k_grs; i++) {
        e1 = quad->kt[i].frent;
        quad->kt[i].frent = map_qd[e1];
        e2 = quad->kt[i].scent;
        if (e2 != -1)
            quad->kt[i].scent = map_qd[e2];
    }
    free(map_qd);

    obs = nd_b;
    temp_nd = (point *) malloc(quad->n_grs * sizeof(point));
    map_nd = (int *) malloc(quad->n_grs * sizeof(int));
    k = 0;
    for (i = 0; i < quad->n_grs; i++)
        if (i != obs) {
            getf_find_rogc_todj(quad->knot[i], &temp_nd[k]);
            map_nd[i] = k;
            z = corresp[i];
            corresp[k] = z;
            k++;
        }

    for (i = 0; i < k; i++)
        getf_find_rogc_todj(temp_nd[i], &quad->knot[i]);
    quad->n_grs = k;
    free(temp_nd);
    for (i = 0; i < quad->e_grs; i++) {
        n1 = quad->elem[i].frvrt;
        quad->elem[i].frvrt = map_nd[n1];
        n2 = quad->elem[i].scvrt;
        quad->elem[i].scvrt = map_nd[n2];
        n3 = quad->elem[i].thvrt;
        quad->elem[i].thvrt = map_nd[n3];
        n4 = quad->elem[i].ftvrt;
        quad->elem[i].ftvrt = map_nd[n4];
    }
    for (i = 0; i < quad->k_grs; i++) {
        n1 = quad->kt[i].frvrt;
        quad->kt[i].frvrt = map_nd[n1];
        n2 = quad->kt[i].scvrt;
        quad->kt[i].scvrt = map_nd[n2];
    }
    free(map_nd);

    inv_ed[n_inv] = new_e_a;
    inv_ed[n_inv + 1] = new_e_b;
    n_inv = n_inv + 2;
    omega = (point *) malloc(max_nb_smooth * sizeof(point));
    type = (int *) malloc(max_nb_smooth * sizeof(int));
    underl = (int *) malloc(max_nb_smooth * sizeof(int));
    for (i = 0; i < n_inv; i++) {
        n1 = quad->kt[inv_ed[i]].frvrt;
        im_n1 = corresp[n1];
        n2 = quad->kt[inv_ed[i]].scvrt;
        im_n2 = corresp[n2];
        nb_stat = folc_gene_kost(im_n1, im_n2, level, nb_inter, GR, msh, omega, type, underl, max_nb_smooth);
        with_conv_davf(msh, GR, omega, type, underl, nb_stat, &fc[inv_ed[i]]);
        if (fc[inv_ed[i]].st_grs >= max_nb_smooth) {
            fprintf(tmpout, "4-Max number of stations [%d]\n", max_nb_smooth);
            exit(0);
        }
        zekl_redu_govr(msh, &fc[inv_ed[i]]);
        nesc_impr_burh(msh, &fc[inv_ed[i]]);
    }
    free(omega);
    free(type);
    free(underl);
    free(inv_ed);

}



void foqn_best_vucp(double t, point * sam, int n_sam, point * X)
{
    int i;
    double *len;
    if (n_sam == 1)
        getf_find_rogc_todj(sam[0], X);
    else {
        len = (double *) malloc(n_sam * sizeof(double));
        len[0] = 0.0;
        for (i = 1; i < n_sam; i++)
            len[i] = len[i - 1] + wodt_dist_gilq(sam[i - 1], sam[i]);
        vefm_eval_bilt(t, sam, n_sam, len, X);
        free(len);
    }
}


void zuvf_best_qitg(double t, float_curve F, point * X)
{
    int n_sam;
    n_sam = F.st_grs;
    foqn_best_vucp(t, F.stn, n_sam, X);
}



void jecr_find_zemc_wozf(fajor_sion3D quad, int nd, int z, float_curve * f, int *e, int *ort)
{
    int ed[4], i, n1, n2, E, N, k;
    double D1, D2;
    point S, T;
    ed[0] = quad.elem[z].frkt;
    ed[1] = quad.elem[z].sckt;
    ed[2] = quad.elem[z].trkt;
    ed[3] = quad.elem[z].ftkt;
    k = 0;
    for (i = 0; i < 4; i++) {
        E = ed[i];
        n1 = quad.kt[E].frvrt;
        n2 = quad.kt[E].scvrt;
        if ((n1 == nd) || (n2 == nd)) {
            N = f[E].st_grs;
            getf_find_rogc_todj(f[E].stn[0], &S);
            getf_find_rogc_todj(f[E].stn[N - 1], &T);
            D1 = wodt_dist_gilq(quad.knot[z], S);
            D2 = wodt_dist_gilq(quad.knot[z], T);
            if (D1 < D2)
                ort[k] = +1;
            else
                ort[k] = -1;
            e[k] = E;
            k++;
        }
    }
    if (k != 2) {
        fprintf(tmpout, "There should be two incident float curves\n");
        exit(0);
    }
}


int moqr_reco_lecn(fajor_sion3D quad, int F, float_curve * fc, int *nd_a, int *nd_b)
{
    int nd[4], op, e[2], ort[2], i, res;
    double *alpha, lambda = 0.2, eps = 0.9, diff1, diff2;
    point C, D;
    nd[0] = quad.elem[F].frvrt;
    nd[1] = quad.elem[F].scvrt;
    nd[2] = quad.elem[F].thvrt;
    nd[3] = quad.elem[F].ftvrt;
    alpha = (double *) malloc(4 * sizeof(double));
    for (i = 0; i < 4; i++) {
        jecr_find_zemc_wozf(quad, nd[i], F, fc, e, ort);
        if (ort[0] == +1)
            zuvf_best_qitg(lambda, fc[e[0]], &C);
        if (ort[0] == -1)
            zuvf_best_qitg(1.0 - lambda, fc[e[0]], &C);
        if (ort[1] == +1)
            zuvf_best_qitg(lambda, fc[e[1]], &D);
        if (ort[1] == -1)
            zuvf_best_qitg(1.0 - lambda, fc[e[1]], &D);
        alpha[i] = bokv_angl_qufn(quad.knot[nd[i]], C, D);
    }
    res = 0;
    for (i = 0; i < 2; i++) {
        if (i == 0)
            op = 2;
        if (i == 1)
            op = 3;
        diff1 = fabs(alpha[i] - MY_PI);
        diff2 = fabs(alpha[op] - MY_PI);
        if ((diff1 < eps) && (diff2 < eps)) {
            res = 1;
            *nd_a = nd[i];
            *nd_b = nd[op];
            break;
        }
    }
    free(alpha);
    return res;
}


void tuqk_trav_nesv(fajor_sion3D * quad, manif_tl msh, prat_main GR, int *corresp, float_curve * fc, int *nb_smooth, int max_nb_smooth)
{
    int nb_inter = 5, i, F, ts, nd_a, nd_b, lev = 0, id, M;
    M = *nb_smooth;
    for (i = 0; i < quad->e_grs; i++) {
        F = i;
        ts = moqr_reco_lecn(*quad, F, fc, &nd_a, &nd_b);
        if (ts == 1) {
            id = cahq_idea_fewg(*quad, corresp, GR, msh, lev, nd_a, nd_b);





            jeqg_skit_leqb(quad, msh, GR, F, nd_a, nd_b, id, corresp, fc, nb_inter, max_nb_smooth);
            M = M - 2;
        }
    }
    *nb_smooth = M;
}
