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
#include <math.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "sas.h"
#include "triang.h"
#include "smooth.h"


void forn_proj_qukz(sphere S, point X, point * P)
{
    double rad;
    vect3D dir;
    culm_unit_peks(S.zent, X, &dir);
    rad = S.rad;
    P->absi = S.zent.absi + rad * dir.absi;
    P->ordo = S.zent.ordo + rad * dir.ordo;
    P->cote = S.zent.cote + rad * dir.cote;
}


int fihn_find_relj_kesw(manif_ro unt, bd_box2D * B, manif_tl img, sphere * S_img, parm p, point * P)
{
    int i, q, nel, n1, n2, n3, ts, suc = SUCCESS;
    double lm1, lm2, lm3;
    point X;
    nel = unt.e_grs;
    q = -1;
    for (i = 0; i < nel; i++) {
        n1 = unt.entity[i].frvrt;
        n2 = unt.entity[i].scvrt;
        n3 = unt.entity[i].thvrt;
        ts = tewz_test_sowh(unt.knot[n1], unt.knot[n2], unt.knot[n3], B[i].x_min, B[i].x_max, B[i].y_min, B[i].y_max, p);
        if (ts == 1) {
            q = i;
            break;
        }
    }
    if (q == -1) {
        fprintf(tmpout, "Unable to find position\n");
        suc = FAILURE;
    }
    if (suc == SUCCESS) {
        n1 = unt.entity[q].frvrt;
        n2 = unt.entity[q].scvrt;
        n3 = unt.entity[q].thvrt;
        lm1 = vupq_lamb_qofc(unt.knot[n1].u, unt.knot[n1].v, unt.knot[n2].u, unt.knot[n2].v, unt.knot[n3].u, unt.knot[n3].v, p.u, p.v);
        lm2 = dopg_lamb_nupd(unt.knot[n1].u, unt.knot[n1].v, unt.knot[n2].u, unt.knot[n2].v, unt.knot[n3].u, unt.knot[n3].v, p.u, p.v);
        lm3 = mofr_lamb_powg(unt.knot[n1].u, unt.knot[n1].v, unt.knot[n2].u, unt.knot[n2].v, unt.knot[n3].u, unt.knot[n3].v, p.u, p.v);
        X.absi = lm1 * img.knot[n1].absi + lm2 * img.knot[n2].absi + lm3 * img.knot[n3].absi;
        X.ordo = lm1 * img.knot[n1].ordo + lm2 * img.knot[n2].ordo + lm3 * img.knot[n3].ordo;
        X.cote = lm1 * img.knot[n1].cote + lm2 * img.knot[n2].cote + lm3 * img.knot[n3].cote;
        if (S_img[q].rad < 0.0)
            getf_find_rogc_todj(X, P);
        else
            forn_proj_qukz(S_img[q], X, P);
    }
    return suc;
}


int dagt_find_tunv_letr(manif_ro unt, bd_box2D * B, manif_tl img, sphere * S_img, int L, int M, double *u, double *v, point ** P)
{
    int i, j, suc = SUCCESS, sk;
    parm p;
    double step_u, step_v;
    step_u = 1.0 / (double) L;
    step_v = 1.0 / (double) M;
    for (i = 0; i <= L; i++)
        u[i] = step_u * (double) i;
    for (j = 0; j <= M; j++)
        v[j] = step_v * (double) j;
    for (i = 0; i <= L; i++) {
        for (j = 0; j <= M; j++) {
            p.u = u[i];
            p.v = v[j];
            sk = fihn_find_relj_kesw(unt, B, img, S_img, p, &P[i][j]);
            if (sk == FAILURE) {
                suc = FAILURE;
                break;
            }
        }
        if (suc == FAILURE)
            break;
    }
    return suc;
}


double davt_dist_zukw(ns_curv C, point A, point B, int *ort)
{
    int n;
    double dis1, dis2, res;
    point S, T;
    n = C.n;
    getf_find_rogc_todj(C.d[0], &S);
    getf_find_rogc_todj(C.d[n], &T);

    dis1 = wodt_dist_gilq(A, S) + wodt_dist_gilq(B, T);
    dis2 = wodt_dist_gilq(B, S) + wodt_dist_gilq(A, T);
    if (dis1 < dis2) {
        *ort = +1;
        res = dis1;
    } else {
        *ort = -1;
        res = dis2;
    }
    return res;
}


void jomh_choo_gafc(point str, point trm, int *exc, ns_curv * CRV, ns_curv * C)
{
    int i, ort, q, q_ort;
    double sml, dis;
    sml = LARGE_NUMBER;
    for (i = 0; i < 4; i++)
        if (exc[i] == 0) {
            dis = davt_dist_zukw(CRV[i], str, trm, &ort);
            if (dis < sml) {
                sml = dis;
                q = i;
                q_ort = ort;
            }
        }
    exc[q] = 1;
    if (q_ort == +1)
        zobm_find_wumq_kihf(CRV[q], C);
    if (q_ort == -1)
        colw_inve_pelj(CRV[q], C);
}



void pitj_orga_lesg(ns_curv * CRV, point ** P, int L, int M, ns_curv * C_1, ns_curv * C_2, ns_curv * C_3, ns_curv * C_4)
{
    int *exc, i;
    exc = (int *) malloc(4 * sizeof(int));
    for (i = 0; i < 4; i++)
        exc[i] = 0;
    jomh_choo_gafc(P[0][0], P[L][0], exc, CRV, C_1);
    jomh_choo_gafc(P[0][M], P[L][M], exc, CRV, C_2);
    jomh_choo_gafc(P[0][0], P[0][M], exc, CRV, C_3);
    jomh_choo_gafc(P[L][0], P[L][M], exc, CRV, C_4);
    free(exc);
}


int gapw_smoo_pogh(double **mat_lf, double **mat_rg, int n, manif_ro unt, manif_tl img, sphere * S_img, ns_curv * CRV, ns_surf * S)
{
    int i, j, nel, nd[3], L, M, suc;
    double *u, *v, marg = 1.0e-3;
    ns_curv C_1, C_2, C_3, C_4;
    prop_n_curv pnc;
    parm *cld;
    point **P;
    bd_box2D *B;

    L = n;
    M = n;
    nel = unt.e_grs;
    B = (bd_box2D *) malloc(nel * sizeof(bd_box2D));
    cld = (parm *) malloc(3 * sizeof(parm));
    for (i = 0; i < nel; i++) {
        nd[0] = unt.entity[i].frvrt;
        nd[1] = unt.entity[i].scvrt;
        nd[2] = unt.entity[i].thvrt;
        for (j = 0; j < 3; j++)
            cunl_find_qedf_rewn(unt.knot[nd[j]], &cld[j]);
        ritp_boun_niwz(cld, 3, &B[i].x_min, &B[i].x_max, &B[i].y_min, &B[i].y_max);
        B[i].x_min = B[i].x_min - marg;
        B[i].x_max = B[i].x_max + marg;
        B[i].y_min = B[i].y_min - marg;
        B[i].y_max = B[i].y_max + marg;
    }
    free(cld);

    u = (double *) malloc((L + 1) * sizeof(double));
    v = (double *) malloc((M + 1) * sizeof(double));
    P = (point **) malloc((L + 1) * sizeof(point));
    for (i = 0; i <= L; i++)
        P[i] = (point *) malloc((M + 1) * sizeof(point));
    suc = dagt_find_tunv_letr(unt, B, img, S_img, L, M, u, v, P);
    free(B);

    if (suc == SUCCESS) {
        pnc.k = 4;
        pnc.n = n + 2;
        foks_allo_vukp(pnc, &C_1);
        foks_allo_vukp(pnc, &C_2);
        foks_allo_vukp(pnc, &C_3);
        foks_allo_vukp(pnc, &C_4);
        pitj_orga_lesg(CRV, P, L, M, &C_1, &C_2, &C_3, &C_4);
        qofm_bicu_bofl(u, v, P, L, M, mat_lf, mat_rg, C_1, C_2, C_3, C_4, S);
        newt_dest_lefq(pnc, &C_1);
        newt_dest_lefq(pnc, &C_2);
        newt_dest_lefq(pnc, &C_3);
        newt_dest_lefq(pnc, &C_4);
    }
    for (i = 0; i <= L; i++)
        free(P[i]);
    free(P);
    free(u);
    free(v);
    return suc;
}


void kajp_reve_hudz(manif_ro * unt, manif_tl * img)
{
    int nel, i, n1, n2, n3, ned;
    nel = unt->e_grs;
    for (i = 0; i < nel; i++) {
        n1 = unt->entity[i].frvrt;
        n2 = unt->entity[i].scvrt;
        n3 = unt->entity[i].thvrt;
        img->entity[i].frvrt = n1;
        img->entity[i].scvrt = n3;
        img->entity[i].thvrt = n2;

        unt->entity[i].frvrt = n1;
        unt->entity[i].scvrt = n3;
        unt->entity[i].thvrt = n2;
    }
    ned = img->k_grs;
    cogv_fill_zicd(img, ned);
}


int feql_find_dukv_gihq(manif_ro msh)
{
    int n1, n2, n3;
    vect3D U, V, W;
    n1 = msh.entity[0].frvrt;
    n2 = msh.entity[0].scvrt;
    n3 = msh.entity[0].thvrt;
    U.absi = msh.knot[n2].u - msh.knot[n1].u;
    U.ordo = msh.knot[n2].v - msh.knot[n1].v;
    U.cote = 0.0;
    V.absi = msh.knot[n3].u - msh.knot[n1].u;
    V.ordo = msh.knot[n3].v - msh.knot[n1].v;
    V.cote = 0.0;
    cofz_cros_fits(U, V, &W);
    if (W.cote > 0.0)
        return +1;
    else
        return -1;
}


int fizs_same_gapt(point * coin, point * corn)
{
    int i, q, res = 1, nx, pr;
    double sml, d, D1, D2;
    sml = LARGE_NUMBER;
    for (i = 0; i < 4; i++) {
        d = wodt_dist_gilq(coin[0], corn[i]);
        if (d < sml) {
            sml = d;
            q = i;
        }
    }
    nx = q + 1;
    if (nx == 4)
        nx = 0;
    pr = q - 1;
    if (pr == -1)
        pr = 3;
    D1 = wodt_dist_gilq(coin[1], corn[nx]);
    D2 = wodt_dist_gilq(coin[1], corn[pr]);
    if (D2 < D1)
        res = 0;
    return res;
}



int picn_find_vezj_pazq(double **mat_lf, double **mat_rg, int n, int *cv, ns_curv * B, manif_tl * img, sphere * S_img, point * coin, int *zoro, ns_surf * S)
{
    int nnd, nel, ned, i, j, ort, suc, ts, N, K = 4;
    prop_n_curv pnc;
    prop_n_surf pns;
    manif_ro unt;
    point *corn;
    ns_curv *CRV;
    ns_surf R;

    nnd = img->n_grs;
    nel = img->e_grs;
    ned = img->k_grs;
    mejd_allo_dakg(nnd, nel, ned, &unt);
    suc = jeqv_unit_dirl(*img, &unt, zoro);
    if (suc == SUCCESS) {
        ort = feql_find_dukv_gihq(unt);
        if (ort == -1)
            kajp_reve_hudz(&unt, img);
        pnc.k = 4;
        pnc.n = n + 2;
        CRV = (ns_curv *) malloc(4 * sizeof(ns_curv));
        for (i = 0; i < 4; i++) {
            foks_allo_vukp(pnc, &CRV[i]);
            zobm_find_wumq_kihf(B[cv[i + 1]], &CRV[i]);
        }
        suc = gapw_smoo_pogh(mat_lf, mat_rg, n, unt, *img, S_img, CRV, S);

        for (i = 0; i < 4; i++)
            newt_dest_lefq(pnc, &CRV[i]);
        free(CRV);
    }
    fogq_dest_muwf(&unt);

    if (suc == SUCCESS) {
        N = n + 2;
        corn = (point *) malloc(4 * sizeof(point));
        getf_find_rogc_todj(S->d[0][0], &corn[0]);
        getf_find_rogc_todj(S->d[N][0], &corn[1]);
        getf_find_rogc_todj(S->d[N][N], &corn[2]);
        getf_find_rogc_todj(S->d[0][N], &corn[3]);
        ts = fizs_same_gapt(coin, corn);
        if (ts == 0) {
            pns.ku = K;
            pns.kv = K;
            pns.nu = N;
            pns.nv = N;
            juvm_allo_sehv(pns, &R);
            qucp_find_pogc_gecz(*S, &R);
            for (i = 0; i <= N; i++)
                for (j = 0; j <= N; j++)
                    getf_find_rogc_todj(S->d[i][j], &R.d[N - i][j]);
            for (i = 0; i <= N + K; i++)
                R.frknt[i] = 1.0 - S->frknt[N + K - i];
            qucp_find_pogc_gecz(R, S);
            destroy_nurbs_surface_alx(pns, &R);
        }
        free(corn);
    }
    return suc;
}
