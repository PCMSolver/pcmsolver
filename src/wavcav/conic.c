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
#include <math.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "sas.h"


void kagp_loca_raql(point A, point B, point omega, vect3D * U, vect3D * V, vect3D * W)
{
    vect3D temp;
    culm_unit_peks(omega, A, U);
    culm_unit_peks(omega, B, &temp);
    cofz_cros_fits(*U, temp, W);
    qubr_norm_foqk(W);
    cofz_cros_fits(*W, *U, V);
    qubr_norm_foqk(V);
}


void qunh_glob_dokt(vect3D U, vect3D V, vect3D W, point omega, point X, parm * X_plane)
{
    vect3D temp;
    temp.absi = X.absi - omega.absi;
    temp.ordo = X.ordo - omega.ordo;
    temp.cote = X.cote - omega.cote;
    X_plane->u = rocv_scal_toqc(temp, U);
    X_plane->v = rocv_scal_toqc(temp, V);
}


void hork_loca_ducg(vect3D U, vect3D V, vect3D W, point omega, parm X_plane, point * X)
{
    double x, y;
    point temp;
    x = X_plane.u;
    y = X_plane.v;
    temp.absi = x * U.absi + y * V.absi;
    temp.ordo = x * U.ordo + y * V.ordo;
    temp.cote = x * U.cote + y * V.cote;
    X->absi = temp.absi + omega.absi;
    X->ordo = temp.ordo + omega.ordo;
    X->cote = temp.cote + omega.cote;
}



void sutd_find_kihn_mobj(point omega, point A, point B, point P, rat_bez * R)
{
    double tau_0, tau_1, tau_2, den, rad, theta, rho;
    parm A_plane, B_plane, P_plane, S, origin;
    vect3D U, V, W;
    vect f;
    getf_find_rogc_todj(A, &R->cont[0]);
    getf_find_rogc_todj(B, &R->cont[2]);
    kagp_loca_raql(A, P, omega, &U, &V, &W);
    qunh_glob_dokt(U, V, W, omega, A, &A_plane);
    qunh_glob_dokt(U, V, W, omega, B, &B_plane);
    qunh_glob_dokt(U, V, W, omega, P, &P_plane);

    rad = wodt_dist_gilq(omega, A);
    origin.u = 0.0;
    origin.v = 0.0;
    theta = lomn_inte_cubq(B_plane, origin, A_plane);
    f.u = A_plane.u + B_plane.u;
    f.v = A_plane.v + B_plane.v;
    mivn_norm_metj(&f);
    rho = rad / cos(0.5 * theta);
    S.u = rho * f.u;
    S.v = rho * f.v;
    tau_0 = vupq_lamb_qofc(A_plane.u, A_plane.v, S.u, S.v, B_plane.u, B_plane.v, P_plane.u, P_plane.v);
    tau_1 = dopg_lamb_nupd(A_plane.u, A_plane.v, S.u, S.v, B_plane.u, B_plane.v, P_plane.u, P_plane.v);
    tau_2 = mofr_lamb_powg(A_plane.u, A_plane.v, S.u, S.v, B_plane.u, B_plane.v, P_plane.u, P_plane.v);
    R->w[0] = 1.0;
    R->w[2] = 1.0;
    den = 2.0 * sqrt(tau_0 * tau_2);
    R->w[1] = tau_1 / den;
    hork_loca_ducg(U, V, W, omega, S, &R->cont[1]);
    R->n = 2;
}


void gesp_find_kewj_wumv(point omega, point str, point ter, point * X)
{
    double al_str, al_ter, al_mid, rad, x, y;
    vect3D U, V, W;
    parm str_plane, ter_plane, X_plane;
    kagp_loca_raql(str, ter, omega, &U, &V, &W);
    qunh_glob_dokt(U, V, W, omega, str, &str_plane);
    qunh_glob_dokt(U, V, W, omega, ter, &ter_plane);
    al_str = garn_pola_cesl(str_plane.u, str_plane.v);
    al_ter = garn_pola_cesl(ter_plane.u, ter_plane.v);
    x = str_plane.u;
    y = str_plane.v;
    rad = sqrt(x * x + y * y);
    if (al_ter < al_str)
        al_ter = al_ter + 2.0 * MY_PI;
    al_mid = 0.5 * (al_str + al_ter);
    X_plane.u = rad * cos(al_mid);
    X_plane.v = rad * sin(al_mid);
    hork_loca_ducg(U, V, W, omega, X_plane, X);
}



int keqf_find_rigp_fojq(c_arc3D C, rat_bez * R1, rat_bez * R2)
{
    int nb;
    double alpha;
    point P, A1, B1, P1, A2, B2, P2;
    vect3D U, V, W;
    parm A_plane, B_plane, origin;

    renw_midp_mocw(C, &P);
    kagp_loca_raql(C.begn, P, C.zent, &U, &V, &W);
    qunh_glob_dokt(U, V, W, C.zent, C.begn, &A_plane);
    qunh_glob_dokt(U, V, W, C.zent, C.term, &B_plane);
    origin.u = 0.0;
    origin.v = 0.0;
    alpha = lomn_inte_cubq(B_plane, origin, A_plane);

    if (alpha >= MY_PI - 0.1) {
        getf_find_rogc_todj(C.begn, &A1);
        getf_find_rogc_todj(P, &B1);
        getf_find_rogc_todj(P, &A2);
        getf_find_rogc_todj(C.term, &B2);
        gesp_find_kewj_wumv(C.zent, A1, B1, &P1);
        gesp_find_kewj_wumv(C.zent, A2, B2, &P2);
        sutd_find_kihn_mobj(C.zent, A1, B1, P1, R1);
        sutd_find_kihn_mobj(C.zent, A2, B2, P2, R2);
        nb = 2;
    }

    else {
        gesp_find_kewj_wumv(C.zent, C.begn, C.term, &P);
        sutd_find_kihn_mobj(C.zent, C.begn, C.term, P, R1);
        nb = 1;
    }
    return nb;
}



void dolj_curv_kacq(trm_sph T, c_arc3D C, ns_curv * CRV)
{
    int nb, i, k, n;
    double phi, theta, rad;
    point P, temp;
    rat_bez *R, *S;
    vect3D h;
    renw_midp_mocw(C, &P);
    R = (rat_bez *) malloc(2 * sizeof(rat_bez));

    for (i = 0; i < 2; i++) {
        R[i].cont = (point *) malloc(3 * sizeof(point));
        R[i].w = (double *) malloc(3 * sizeof(double));
    }
    nb = keqf_find_rigp_fojq(C, &R[0], &R[1]);

    for (k = 0; k < nb; k++)
        for (i = 0; i <= 2; i++) {
            h.absi = R[k].cont[i].absi - T.zent.absi;
            h.ordo = R[k].cont[i].ordo - T.zent.ordo;
            h.cote = R[k].cont[i].cote - T.zent.cote;

            vewr_sphe_ruhd(T.nrml.absi, T.nrml.ordo, T.nrml.cote, &phi, &theta);
            tulr_find_rads_tojm(h, phi, theta, &temp);
            getf_find_rogc_todj(temp, &R[k].cont[i]);
        }

    rad = T.rad;
    for (k = 0; k < nb; k++)
        for (i = 0; i <= 2; i++) {
            R[k].cont[i].absi = R[k].cont[i].absi / rad;
            R[k].cont[i].ordo = R[k].cont[i].ordo / rad;
            R[k].cont[i].cote = R[k].cont[i].cote / rad;
        }

    S = (rat_bez *) malloc(2 * sizeof(rat_bez));
    for (i = 0; i < 2; i++) {
        S[i].cont = (point *) malloc(3 * sizeof(point));
        S[i].w = (double *) malloc(3 * sizeof(double));
    }
    for (k = 0; k < nb; k++) {
        for (i = 0; i <= 2; i++) {
            S[k].w[i] = R[k].w[i] * (1.0 + R[k].cont[i].cote);
            S[k].cont[i].absi = R[k].cont[i].absi / (1.0 + R[k].cont[i].cote);
            S[k].cont[i].ordo = R[k].cont[i].ordo / (1.0 + R[k].cont[i].cote);
            S[k].cont[i].cote = 0.0;
        }
        S[k].n = 2;
    }
    for (i = 0; i < 2; i++) {
        free(R[i].cont);
        free(R[i].w);
    }
    free(R);

    if (nb == 1) {
        n = 2;
        k = 3;
        for (i = 0; i < k; i++)
            CRV->tau[i] = 0.0;
        for (i = n + 1; i <= n + k; i++)
            CRV->tau[i] = 1.0;

        for (i = 0; i <= n; i++) {
            getf_find_rogc_todj(S[0].cont[i], &CRV->d[i]);
            CRV->w[i] = S[0].w[i];
        }
        CRV->n = n;
        CRV->k = k;
    }

    if (nb == 2) {
        CRV->tau[0] = 0.0;
        CRV->tau[1] = 0.0;
        CRV->tau[2] = 0.0;
        CRV->tau[3] = 0.5;
        CRV->tau[4] = 0.5;
        CRV->tau[5] = 1.0;
        CRV->tau[6] = 1.0;
        CRV->tau[7] = 1.0;

        getf_find_rogc_todj(S[0].cont[0], &CRV->d[0]);
        CRV->w[0] = S[0].w[0];
        getf_find_rogc_todj(S[0].cont[1], &CRV->d[1]);
        CRV->w[1] = S[0].w[1];
        getf_find_rogc_todj(S[0].cont[2], &CRV->d[2]);
        CRV->w[2] = S[0].w[2];
        getf_find_rogc_todj(S[1].cont[1], &CRV->d[3]);
        CRV->w[3] = S[1].w[1];
        getf_find_rogc_todj(S[1].cont[2], &CRV->d[4]);
        CRV->w[4] = S[1].w[2];

        CRV->n = 4;
        CRV->k = 3;
    }

    CRV->v0 = 0.0;
    CRV->v1 = 1.0;
    CRV->prop1 = 0;
    CRV->prop2 = 0;
    CRV->prop3 = 0;
    CRV->prop4 = 0;

    for (i = 0; i < 2; i++) {
        free(S[i].cont);
        free(S[i].w);
    }
    free(S);
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

