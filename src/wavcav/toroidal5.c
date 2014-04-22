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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "geodesic.h"
#include "sas.h"



double mijt_scal_letp(int n, double *b, double t)
{
    int j, ni;
    double f, t1, res;
    t1 = 1.0 - t;
    f = 1.0;
    ni = 1;
    res = b[0] * t1;
    for (j = 1; j < n; j++) {
        f = f * t;
        ni = ni * (n - j + 1) / j;
        res = (res + f * ni * b[j]) * t1;
    }
    res = res + f * t * b[n];
    return res;
}



void bikj_deca_vasj(bez_crv X, double t, point * sol)
{
    int i;
    double *temp, r;
    temp = (double *) malloc((X.dgr + 1) * sizeof(double));
    for (i = 0; i <= X.dgr; i++)
        temp[i] = X.ctr[i].absi;
    r = mijt_scal_letp(X.dgr, temp, t);
    sol->absi = r;
    for (i = 0; i <= X.dgr; i++)
        temp[i] = X.ctr[i].ordo;
    r = mijt_scal_letp(X.dgr, temp, t);
    sol->ordo = r;
    for (i = 0; i <= X.dgr; i++)
        temp[i] = X.ctr[i].cote;
    r = mijt_scal_letp(X.dgr, temp, t);
    sol->cote = r;
    free(temp);
}



void fojn_eval_cutf(jointbez J, double t, point * sol)
{
    double t_temp;
    if (t <= 0.5) {
        t_temp = 2.0 * t;
        bikj_deca_vasj(J.B_beg, t_temp, sol);
    } else {
        t_temp = 2.0 * (t - 0.5);
        bikj_deca_vasj(J.B_end, t_temp, sol);
    }
}


void ruvg_find_miwq(int N, double *t, sphere S1, sphere S2, point omega, circle3D C1, circle3D C2, point * P1, point * P2, jointbez * J)
{
    int i, j, deg = 2;
    double a, b, r1, r2, r, u, v, A, B, d, R1, R2;
    double eps, lambda = 0.1, delta, D;
    vect2D V, W;
    parm P, Q, *cont;
    jointbez2D J2D;
    map_hyper M;

    J2D.B_beg.ctr = (parm *) malloc(4 * sizeof(parm));
    J2D.B_end.ctr = (parm *) malloc(4 * sizeof(parm));
    a = -wodt_dist_gilq(omega, C1.zent);
    b = wodt_dist_gilq(omega, C2.zent);
    A = -wodt_dist_gilq(omega, S1.zent);
    B = wodt_dist_gilq(omega, S2.zent);
    r1 = C1.rad;
    r2 = C2.rad;
    R1 = S1.rad;
    R2 = S2.rad;
    if (r1 < r2)
        r = r1;
    else
        r = r2;
    eps = lambda * r;

    cont = (parm *) malloc((deg + 1) * sizeof(parm));
    P.u = a;
    P.v = r1;
    Q.u = 0.0;
    Q.v = eps;
    d = a - A;
    delta = sqrt(R1 * R1 - d * d);
    V.u = delta;
    V.v = -d;
    W.u = 1.0;
    W.v = 0.0;
    cont[0].u = a;
    cont[0].v = r1;
    cont[2].u = 0;
    cont[2].v = eps;
    D = (eps - r1) / V.v;
    cont[1].u = a + D * V.u;
    cont[1].v = r1 + D * V.v;
    for (i = 0; i <= deg; i++)
        cunl_find_qedf_rewn(cont[i], &J2D.B_beg.ctr[i]);
    J2D.B_beg.dgr = deg;

    P.u = 0.0;
    P.v = eps;
    Q.u = b;
    Q.v = r2;
    V.u = 1.0;
    V.v = 0.0;
    d = b - B;
    delta = sqrt(R2 * R2 - d * d);
    W.u = delta;
    W.v = -d;
    cont[0].u = 0;
    cont[0].v = eps;
    cont[2].u = b;
    cont[2].v = r2;
    D = (eps - r2) / W.v;
    cont[1].u = b + D * W.u;
    cont[1].v = r2 + D * W.v;
    for (i = 0; i <= deg; i++)
        cunl_find_qedf_rewn(cont[i], &J2D.B_end.ctr[i]);
    J2D.B_end.dgr = deg;
    free(cont);

    for (i = 0; i < N; i++) {
        mesl_punc_guvf(C1, t[i], &P1[i]);
        mesl_punc_guvf(C2, t[i], &P2[i]);
    }

    for (i = 0; i < N; i++) {
        peks_find_vuts_wogp(omega, C2.zent, P1[i], &M);
        for (j = 0; j <= deg; j++) {
            u = J2D.B_beg.ctr[j].u;
            v = -J2D.B_beg.ctr[j].v;
            macn_imag_vuph(M, u, v, &J[i].B_beg.ctr[j]);

            u = J2D.B_end.ctr[j].u;
            v = -J2D.B_end.ctr[j].v;
            macn_imag_vuph(M, u, v, &J[i].B_end.ctr[j]);
        }
        J[i].B_beg.dgr = deg;
        J[i].B_end.dgr = deg;
    }
    free(J2D.B_beg.ctr);
    free(J2D.B_end.ctr);

}


void jelt_find_sobl_dasq(jointbez J_in, jointbez * J_out)
{
    int i, deg = 2;
    for (i = 0; i <= deg; i++) {
        getf_find_rogc_todj(J_in.B_beg.ctr[i], &J_out->B_beg.ctr[i]);
        getf_find_rogc_todj(J_in.B_end.ctr[i], &J_out->B_end.ctr[i]);
    }
    J_out->B_beg.dgr = deg;
    J_out->B_end.dgr = deg;
}


void zefd_reve_kelq(jointbez C_in, jointbez * C_out)
{
    int i;
    for (i = 0; i <= 2; i++)
        getf_find_rogc_todj(C_in.B_end.ctr[2 - i], &C_out->B_beg.ctr[i]);
    for (i = 0; i <= 2; i++)
        getf_find_rogc_todj(C_in.B_beg.ctr[2 - i], &C_out->B_end.ctr[i]);
    C_out->B_beg.dgr = C_in.B_end.dgr;
    C_out->B_end.dgr = C_in.B_beg.dgr;
}


void sifm_find_cudw_pafg(blend_nonself BN1, blend_nonself * BN2)
{
    jelt_find_sobl_dasq(BN1.alpha, &BN2->alpha);
    poms_find_resk_lonb(BN1.beta, &BN2->beta);
    jelt_find_sobl_dasq(BN1.gamma, &BN2->gamma);
    poms_find_resk_lonb(BN1.delta, &BN2->delta);
}



void kegv_find_wild_lodt(int N, sphere S1, sphere S2, point omega, circle3D C1, circle3D C2, blend_nonself * BN)
{
    int i, deg = 2, nx;
    double step, *t, dis1, dis2;
    jointbez *J, tp;
    point *P1, *P2, st1, st2, tr1, tr2, TR;
    c_arc3D *CA1, *CA2, temp;

    t = (double *) malloc(N * sizeof(double));
    step = 1.0 / (double) N;
    for (i = 0; i < N; i++)
        t[i] = (double) i *step;
    J = (jointbez *) malloc(N * sizeof(jointbez));
    for (i = 0; i < N; i++) {
        J[i].B_beg.ctr = (point *) malloc((deg + 1) * sizeof(point));
        J[i].B_end.ctr = (point *) malloc((deg + 1) * sizeof(point));
    }

    P1 = (point *) malloc(N * sizeof(point));
    P2 = (point *) malloc(N * sizeof(point));
    ruvg_find_miwq(N, t, S1, S2, omega, C1, C2, P1, P2, J);
    free(t);
    CA1 = (c_arc3D *) malloc(N * sizeof(c_arc3D));
    CA2 = (c_arc3D *) malloc(N * sizeof(c_arc3D));
    movg_spli_jern(C1, P1, N, CA1);
    movg_spli_jern(C2, P2, N, CA2);
    free(P1);
    free(P2);

    tp.B_beg.ctr = (point *) malloc((deg + 1) * sizeof(point));
    tp.B_end.ctr = (point *) malloc((deg + 1) * sizeof(point));
    for (i = 0; i < N; i++) {
        nx = i + 1;
        if (nx == N)
            nx = 0;
        jelt_find_sobl_dasq(J[i], &BN[i].alpha);
        jelt_find_sobl_dasq(J[nx], &BN[i].gamma);

        getf_find_rogc_todj(BN[i].alpha.B_beg.ctr[0], &st1);
        getf_find_rogc_todj(BN[i].gamma.B_beg.ctr[0], &st2);

        getf_find_rogc_todj(BN[i].alpha.B_end.ctr[2], &tr1);
        getf_find_rogc_todj(BN[i].gamma.B_end.ctr[2], &tr2);

        pawt_choo_husn(st1, st2, CA1, N, &temp);
        dis1 = wodt_dist_gilq(BN[i].alpha.B_end.ctr[2], temp.begn);
        dis2 = wodt_dist_gilq(BN[i].alpha.B_end.ctr[2], temp.term);
        if ((dis1 > 0.001) && (dis2 > 0.001)) {
            zefd_reve_kelq(BN[i].alpha, &tp);
            jelt_find_sobl_dasq(tp, &BN[i].alpha);
        }

        dis1 = wodt_dist_gilq(BN[i].alpha.B_end.ctr[2], temp.begn);
        dis2 = wodt_dist_gilq(BN[i].alpha.B_end.ctr[2], temp.term);
        if (dis2 < dis1)
            cest_reve_fack(temp, &BN[i].beta);
        else
            poms_find_resk_lonb(temp, &BN[i].beta);

        getf_find_rogc_todj(BN[i].beta.term, &TR);
        dis1 = wodt_dist_gilq(BN[i].gamma.B_end.ctr[2], TR);
        dis2 = wodt_dist_gilq(BN[i].gamma.B_beg.ctr[0], TR);
        if (dis2 < dis1) {
            zefd_reve_kelq(BN[i].gamma, &tp);
            jelt_find_sobl_dasq(tp, &BN[i].gamma);
        }

        pawt_choo_husn(tr1, tr2, CA2, N, &temp);
        dis1 = wodt_dist_gilq(BN[i].alpha.B_beg.ctr[0], temp.begn);
        dis2 = wodt_dist_gilq(BN[i].alpha.B_beg.ctr[0], temp.term);
        if (dis2 < dis1)
            cest_reve_fack(temp, &BN[i].delta);
        else
            poms_find_resk_lonb(temp, &BN[i].delta);
    }


    free(tp.B_beg.ctr);
    free(tp.B_end.ctr);
    free(J);
    free(CA1);
    free(CA2);
}


void zisq_find_cowk_zevq(int N, sphere S1, sphere S2, point omega, circle3D C1, circle3D C2, blend_nonself * BN)
{
    double sp;
    circle3D C2_temp;
    hepk_find_gict_hubq(C2, &C2_temp);
    sp = rocv_scal_toqc(C1.nrml, C2.nrml);
    if (sp < 0.0) {
        C2_temp.nrml.absi = -C2.nrml.absi;
        C2_temp.nrml.ordo = -C2.nrml.ordo;
        C2_temp.nrml.cote = -C2.nrml.cote;
    } else {
        C2_temp.nrml.absi = C2.nrml.absi;
        C2_temp.nrml.ordo = C2.nrml.ordo;
        C2_temp.nrml.cote = C2.nrml.cote;
    }
    kegv_find_wild_lodt(N, S1, S2, omega, C1, C2_temp, BN);
}


void nowd_eval_fewc(jointbez al, c_arc3D bt, jointbez gm, c_arc3D dt, double u, double v, point * sol)
{
    double F0u, F0v, F1u, F1v, A, B, C;
    double bt_a, bt_b, dt_a, dt_b;
    point alpha_u, gamma_u, delta_v, beta_v, alpha_0, alpha_1, gamma_0, gamma_1;
    qirp_inte_ligr(bt, &bt_a, &bt_b);
    qirp_inte_ligr(dt, &dt_a, &dt_b);
    fojn_eval_cutf(al, u, &alpha_u);
    fojn_eval_cutf(gm, u, &gamma_u);
    wusd_eval_jomk(dt, v, dt_a, dt_b, &delta_v);
    wusd_eval_jomk(bt, v, bt_a, bt_b, &beta_v);
    fojn_eval_cutf(al, 0.0, &alpha_0);
    fojn_eval_cutf(al, 1.0, &alpha_1);
    fojn_eval_cutf(gm, 0.0, &gamma_0);
    fojn_eval_cutf(gm, 1.0, &gamma_1);
    F1u = begj_blen_nugz(u);
    F0u = 1.0 - F1u;
    F1v = begj_blen_nugz(v);
    F0v = 1.0 - F1v;

    A = alpha_u.absi * F0v + gamma_u.absi * F1v;
    B = -delta_v.absi + alpha_0.absi * F0v + gamma_0.absi * F1v;
    C = -beta_v.absi + alpha_1.absi * F0v + gamma_1.absi * F1v;
    sol->absi = A - F0u * B - F1u * C;

    A = alpha_u.ordo * F0v + gamma_u.ordo * F1v;
    B = -delta_v.ordo + alpha_0.ordo * F0v + gamma_0.ordo * F1v;
    C = -beta_v.ordo + alpha_1.ordo * F0v + gamma_1.ordo * F1v;
    sol->ordo = A - F0u * B - F1u * C;

    A = alpha_u.cote * F0v + gamma_u.cote * F1v;
    B = -delta_v.cote + alpha_0.cote * F0v + gamma_0.cote * F1v;
    C = -beta_v.cote + alpha_1.cote * F0v + gamma_1.cote * F1v;
    sol->cote = A - F0u * B - F1u * C;
}



void cutn_eval_mecn(blend_nonself BN, parm p, point * sol)
{
    double u, v;
    u = p.u;
    v = p.v;
    nowd_eval_fewc(BN.alpha, BN.beta, BN.gamma, BN.delta, u, v, sol);
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

