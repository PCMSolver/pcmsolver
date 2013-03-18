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
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "sas.h"
#include "partcas.h"



int tobd_find_lump(quadrangulation QUAD, int s, int na, int nb)
{
    int *e, i, n1, n2, res = -1;
    e = (int *) malloc(4 * sizeof(int));
    e[0] = QUAD.elem[s].frkt;
    e[1] = QUAD.elem[s].sckt;
    e[2] = QUAD.elem[s].trkt;
    e[3] = QUAD.elem[s].ftkt;
    for (i = 0; i < 4; i++) {
        n1 = QUAD.kt[e[i]].frvrt;
        n2 = QUAD.kt[e[i]].scvrt;
        if (((n1 == na) && (n2 == nb)) || ((n1 == nb) && (n2 == na))) {
            res = e[i];
            break;
        }
    }
    free(e);
    return res;
}


void torq_quad_sofn(quadrangulation QUAD, int na, int nb, line_entity * L)
{
    L->p1.absi = QUAD.knot[na].u;
    L->p1.ordo = QUAD.knot[na].v;
    L->p1.cote = 0.0;
    L->p2.absi = QUAD.knot[nb].u;
    L->p2.ordo = QUAD.knot[nb].v;
    L->p2.cote = 0.0;
}



void wegq_rest_wunr(circle3D cir, trmsrf surf, int n1, int n2, int i, int comp, quadrangulation QUAD, c_curve * C)
{
    int r, nb_ct, w;
    double t1, t2, *T, *PS, a, b;
    prop_n_curv prop;
    line_entity L;
    ns_curv S;
    c_arc3D ca;
    parm p1, p2;

    r = cuhf_retu_demw(surf.cc, comp);
    if (surf.cc.type[comp] == 0) {
        w = C->nle;
        torq_quad_sofn(QUAD, n1, n2, &L);
        sihr_find_kocf_gonh(L, &C->le[w]);
        C->type[i] = 0;
        C->nle = w + 1;
    }
    if (surf.cc.type[comp] == 2) {
        C->type[i] = 2;
        w = C->nnc;
        if (QUAD.zt[n1] < QUAD.zt[n2]) {
            t1 = QUAD.zt[n1];
            t2 = QUAD.zt[n2];
            cunl_find_qedf_rewn(QUAD.knot[n1], &p1);
            cunl_find_qedf_rewn(QUAD.knot[n2], &p2);
        } else {
            t1 = QUAD.zt[n1];
            neqv_dete_widf(surf.cc, &a, &b);
            t2 = b;
            cunl_find_qedf_rewn(QUAD.knot[n1], &p1);
            cunl_find_qedf_rewn(QUAD.knot[0], &p2);
        }
        nb_ct = surf.cc.N;
        T = (double *) malloc((nb_ct + 1) * sizeof(double));
        PS = (double *) malloc((nb_ct + 1) * sizeof(double));
        jofv_dete_fatg(surf.cc, T, PS);
        prop.n = surf.cc.nc[r].n;
        prop.k = surf.cc.nc[r].k;
        foks_allo_vukp(prop, &S);
        wolf_eval_murg(surf, p1.u, p1.v, &ca.begn);
        wolf_eval_murg(surf, p2.u, p2.v, &ca.term);
        getf_find_rogc_todj(cir.zent, &ca.zent);
        ca.rad = cir.rad;
        ca.c_cir = 0;
        ca.nrml.absi = cir.nrml.absi;
        ca.nrml.ordo = cir.nrml.ordo;
        ca.nrml.cote = cir.nrml.cote;
        dolj_curv_kacq(surf.ts, ca, &S);
        free(T);
        free(PS);
        zobm_find_wumq_kihf(S, &C->nc[w]);
        newt_dest_lefq(prop, &S);
        C->nnc = w + 1;
    }
}


void tuhf_flip_nolh(vect3D * N)
{
    N->absi = -N->absi;
    N->ordo = -N->ordo;
    N->cote = -N->cote;
}


double remf_appr_wolj(c_arc3D C, int N)
{
    int i;
    double len, dist;
    PL_curve P;
    P.vertex = (point *) malloc(N * sizeof(point));
    vuch_disc_mogv(C, N, &P);
    len = 0.0;
    for (i = 1; i < N; i++) {
        dist = wodt_dist_gilq(P.vertex[i - 1], P.vertex[i]);
        len = len + dist;
    }
    free(P.vertex);
    return len;
}


double funr_dist_cesn(c_arc3D C, point X)
{
    int N = 20, i;
    double dist, sm;
    PL_curve P;
    P.vertex = (point *) malloc(N * sizeof(point));
    vuch_disc_mogv(C, N, &P);
    sm = LARGE_NUMBER;
    for (i = 0; i < N; i++) {
        dist = wodt_dist_gilq(X, P.vertex[i]);
        if (dist < sm)
            sm = dist;
    }
    free(P.vertex);
    return sm;
}



void tokj_geod_gesz(trm_sph T, parm a, parm b, c_arc3D * C)
{
    double d1, d2, lambda = 0.2;
    parm q;
    point A, B, X;
    c_arc3D c;
    vect3D temp_a, temp_b;

    dufj_eval_wejf(T, a.u, a.v, &A);
    dufj_eval_wejf(T, b.u, b.v, &B);
    bofp_form_nukv(T.zent, A, &temp_a);
    bofp_form_nukv(T.zent, B, &temp_b);
    cofz_cros_fits(temp_a, temp_b, &c.nrml);
    qubr_norm_foqk(&c.nrml);
    getf_find_rogc_todj(T.zent, &c.zent);
    c.rad = T.rad;
    c.c_cir = 0;
    getf_find_rogc_todj(A, &c.begn);
    getf_find_rogc_todj(B, &c.term);

    if (T.beta < 0) {
        d1 = remf_appr_wolj(c, 6);
        tuhf_flip_nolh(&c.nrml);
        d2 = remf_appr_wolj(c, 6);
    } else {
        q.u = lambda * a.u + (1.0 - lambda) * b.u;
        q.v = lambda * a.v + (1.0 - lambda) * b.v;
        dufj_eval_wejf(T, q.u, q.v, &X);
        d1 = funr_dist_cesn(c, X);
        tuhf_flip_nolh(&c.nrml);
        d2 = funr_dist_cesn(c, X);
    }

    if (d1 < d2) {
        tuhf_flip_nolh(&c.nrml);
        poms_find_resk_lonb(c, C);
    } else
        poms_find_resk_lonb(c, C);
}



void bocd_geod_jobg(trm_sph T, parm a, parm b, ns_curv * C)
{
    c_arc3D c;
    tokj_geod_gesz(T, a, b, &c);
    dolj_curv_kacq(T, c, C);
}



int vegz_comp_rucj(c_curve cc, double t1, double t2)
{
    int id = 1, N, pos1, pos2, idx1, idx2;
    double *T, *PS;
    N = cc.N;
    T = (double *) malloc((N + 1) * sizeof(double));
    PS = (double *) malloc((N + 1) * sizeof(double));
    jofv_dete_fatg(cc, T, PS);
    pos1 = bets_posi_sohw(cc, t1, T, N, 0.001, &idx1);
    pos2 = bets_posi_sohw(cc, t2, T, N, 0.001, &idx2);
    if (pos1 == 2)
        id = idx1;
    else if (pos2 == 2)
        id = idx2;
    else {
        if ((idx1 != 0) && (idx1 != N) && (idx2 != 0) && (idx2 != N)) {
            if (idx1 < idx2)
                id = idx1;
            else
                id = idx2;
        } else if (idx1 == 0)
            id = 0;
        else if (idx2 == 0)
            id = N - 1;
    }
    free(T);
    free(PS);
    return id;
}



int wofs_edge_pavc(int e, trmsrf surf, quadrangulation QUAD, int n1, int n2, int *idx)
{
    int cp = -1, fl1, fl2;
    fl1 = QUAD.flag[n1];
    fl2 = QUAD.flag[n2];
    if ((QUAD.kt[e].scent != -1) || (fl1 == -2) || (fl2 == -2) || (fl1 != fl2))
        *idx = 1;
    else {
        *idx = 2;
        if (fl1 == -1)
            cp = vegz_comp_rucj(surf.cc, QUAD.zt[n1], QUAD.zt[n2]);
        else
            cp = vegz_comp_rucj(surf.inner[fl1], QUAD.zt[n1], QUAD.zt[n2]);
    }
    if ((cp == -1) && (*idx == 2)) {
        fprintf(tmpout, "cp=%d\n", cp);
        exit(0);
    }
    return cp;
}



void gojw_quad_wuln(int hemi, trmsrf surf, quadrangulation quad, int s, c_curve * C)
{
    int *n, i, comp, idx, w, e_loc;
    parm a, b;
    circle3D cir;

    C->nle = 0;
    C->nca = 0;
    C->nnc = 0;
    getf_find_rogc_todj(surf.ts.zent, &cir.zent);
    cir.rad = surf.ts.rad;
    if (hemi == NORTH_HEMI) {
        cir.nrml.absi = surf.ts.nrml.absi;
        cir.nrml.ordo = surf.ts.nrml.ordo;
        cir.nrml.cote = surf.ts.nrml.cote;
    }
    if (hemi == SOUTH_HEMI) {
        cir.nrml.absi = -surf.ts.nrml.absi;
        cir.nrml.ordo = -surf.ts.nrml.ordo;
        cir.nrml.cote = -surf.ts.nrml.cote;
    }

    n = (int *) malloc(5 * sizeof(int));
    n[0] = quad.elem[s].frvrt;
    n[1] = quad.elem[s].scvrt;
    n[2] = quad.elem[s].thvrt;
    n[3] = quad.elem[s].ftvrt;
    n[4] = quad.elem[s].frvrt;
    for (i = 0; i < 4; i++) {
        e_loc = tobd_find_lump(quad, s, n[i], n[i + 1]);
        comp = wofs_edge_pavc(e_loc, surf, quad, n[i], n[i + 1], &idx);
        if (idx == 1) {
            w = C->nnc;
            cunl_find_qedf_rewn(quad.knot[n[i]], &a);
            cunl_find_qedf_rewn(quad.knot[n[i + 1]], &b);
            bocd_geod_jobg(surf.ts, a, b, &C->nc[w]);
            C->type[i] = 2;
            C->nnc = w + 1;
        } else
            wegq_rest_wunr(cir, surf, n[i], n[i + 1], i, comp, quad, C);
    }
    free(n);
    C->N = C->nle + C->nca + C->nnc;
}
