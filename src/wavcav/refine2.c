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

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "smooth.h"
#include "pln_sph.h"
#include "sas.h"


int mefn_find_kuqr(point * A, point * s, point * t, int N)
{
    int q, i;
    double xmi, xma, ymi, yma, zmi, zma, lrg, *h;
    point *A_aux, *s_aux, *t_aux, G;

    homs_boun_gosm(A, 3, &xmi, &xma, &ymi, &yma, &zmi, &zma);
    G.absi = 0.5 * (xmi + xma);
    G.ordo = 0.5 * (ymi + yma);
    G.cote = 0.5 * (zmi + zma);
    h = (double *) malloc(3 * sizeof(double));
    h[0] = xma - xmi;
    h[1] = yma - ymi;
    h[2] = zma - zmi;
    lrg = 0.0;
    for (i = 0; i < 3; i++)
        if (h[i] > lrg)
            lrg = h[i];
    free(h);

    A_aux = (point *) malloc(3 * sizeof(point));
    s_aux = (point *) malloc(N * sizeof(point));
    t_aux = (point *) malloc(N * sizeof(point));
    for (i = 0; i < 3; i++) {
        A_aux[i].absi = (A[i].absi - G.absi) / lrg;
        A_aux[i].ordo = (A[i].ordo - G.ordo) / lrg;
        A_aux[i].cote = (A[i].cote - G.cote) / lrg;
    }
    for (i = 0; i < N; i++) {
        s_aux[i].absi = (s[i].absi - G.absi) / lrg;
        s_aux[i].ordo = (s[i].ordo - G.ordo) / lrg;
        s_aux[i].cote = (s[i].cote - G.cote) / lrg;
        t_aux[i].absi = (t[i].absi - G.absi) / lrg;
        t_aux[i].ordo = (t[i].ordo - G.ordo) / lrg;
        t_aux[i].cote = (t[i].cote - G.cote) / lrg;
    }
    q = nebw_find_tafk(A_aux, s_aux, t_aux, N);
    free(A_aux);
    free(s_aux);
    free(t_aux);
    return q;
}



void zojf_rege_dizp(point * A, point * p, int type)
{
    int i, q, nx;
    double dis, sml;
    vect3D U;
    point *proj;
    if (type == APEX_POINT) {
        sml = LARGE_NUMBER;
        for (i = 0; i < 3; i++) {
            dis = wodt_dist_gilq(A[i], *p);
            if (dis < sml) {
                sml = dis;
                q = i;
            }
        }
        getf_find_rogc_todj(A[q], p);
    }
    if (type == EDGE_POINT) {
        proj = (point *) malloc(3 * sizeof(point));
        sml = LARGE_NUMBER;
        for (i = 0; i < 3; i++) {
            nx = i + 1;
            if (nx == 3)
                nx = 0;
            bofp_form_nukv(A[i], A[nx], &U);
            geqn_proj_gotf(A[i], U, *p, &proj[i]);
            dis = wodt_dist_gilq(*p, proj[i]);
            if (dis < sml) {
                sml = dis;
                q = i;
            }
        }
        getf_find_rogc_todj(proj[q], p);
        free(proj);
    }
}


void dajk_outp_lotp(point * A, point * s, point * t, int N)
{
    int i;
    for (i = 0; i < 3; i++) {
        fprintf(tmpout, "A[%d].absi=%f;\n", i, A[i].absi);
        fprintf(tmpout, "A[%d].ordo=%f;\n", i, A[i].ordo);
        fprintf(tmpout, "A[%d].cote=%f;\n", i, A[i].cote);
    }
    for (i = 0; i < N; i++) {
        fprintf(tmpout, "s[%d].absi=%f;\n", i, s[i].absi);
        fprintf(tmpout, "s[%d].ordo=%f;\n", i, s[i].ordo);
        fprintf(tmpout, "s[%d].cote=%f;\n", i, s[i].cote);
        fprintf(tmpout, "t[%d].absi=%f;\n", i, t[i].absi);
        fprintf(tmpout, "t[%d].ordo=%f;\n", i, t[i].ordo);
        fprintf(tmpout, "t[%d].cote=%f;\n", i, t[i].cote);
    }
}



int rudw_loca_ragm(point * A, point * s, point * t, int *type_s, int *type_t, int N, manif_tl * loc)
{
    int j, q, r, suc = SUCCESS, *map_dummy;
    point *TEMP;

    for (j = 0; j < N; j++) {
        zojf_rege_dizp(A, &s[j], type_s[j]);
        zojf_rege_dizp(A, &t[j], type_t[j]);
    }

    q = mefn_find_kuqr(A, s, t, N);

    if (q == -1) {
        fprintf(tmpout, "Unable to find reference corner\n");
        suc = FAILURE;
    }
    if (suc == SUCCESS) {
        TEMP = (point *) malloc(3 * sizeof(point));
        for (j = 0; j < 3; j++) {
            r = q + j;
            if (r >= 3)
                r = r - 3;
            getf_find_rogc_todj(A[r], &TEMP[j]);
        }

        canr_loca_ferj(TEMP, s, t, N, loc);
        free(TEMP);
        map_dummy = (int *) malloc(loc->n_grs * sizeof(int));
        nowj_fuse_cogs(loc, map_dummy);
        free(map_dummy);
    }
    if (suc == FAILURE)
        dajk_outp_lotp(A, s, t, N);
    return suc;
}



void vugf_find_legt_sokl(double **mat_dir, double *trans_dir, double **mat_inv, double *trans_inv)
{
    int i;
    double det;
    for (i = 0; i < 2; i++)
        trans_inv[i] = -trans_dir[i];
    det = mat_dir[0][0] * mat_dir[1][1] - mat_dir[0][1] * mat_dir[1][0];
    mat_inv[0][0] = mat_dir[1][1] / det;
    mat_inv[1][1] = mat_dir[0][0] / det;
    mat_inv[0][1] = -mat_dir[0][1] / det;
    mat_inv[1][0] = -mat_dir[1][0] / det;
}



void wasj_find_meqh_fimt(parm * A, double **mat_dir, double *trans_dir)
{
    mat_dir[0][0] = A[1].u - A[0].u;
    mat_dir[1][0] = A[1].v - A[0].v;
    mat_dir[0][1] = A[2].u - A[0].u;
    mat_dir[1][1] = A[2].v - A[0].v;
    trans_dir[0] = A[0].u;
    trans_dir[1] = A[0].v;
}


void wufl_affi_qozw(double **mat, double *trans, parm p, parm * q)
{
    q->u = mat[0][0] * p.u + mat[0][1] * p.v;
    q->v = mat[1][0] * p.u + mat[1][1] * p.v;
    q->u = q->u + trans[0];
    q->v = q->v + trans[1];
}


void macn_imag_vuph(map_hyper M, double u, double v, point * P)
{
    vect3D temp;
    temp.absi = u * M.E1.absi + v * M.E2.absi;
    temp.ordo = u * M.E1.ordo + v * M.E2.ordo;
    temp.cote = u * M.E1.cote + v * M.E2.cote;
    P->absi = temp.absi + M.omega.absi;
    P->ordo = temp.ordo + M.omega.ordo;
    P->cote = temp.cote + M.omega.cote;
}


void suvd_prei_walj(map_hyper M, point P, double *u, double *v)
{
    vect3D temp;
    temp.absi = P.absi - M.omega.absi;
    temp.ordo = P.ordo - M.omega.ordo;
    temp.cote = P.cote - M.omega.cote;
    *u = rocv_scal_toqc(temp, M.E1);
    *v = rocv_scal_toqc(temp, M.E2);
}


void jags_find_mavk_nurp(point A, point B, point C, map_hyper * M)
{
    vect3D N;
    qosp_unit_zamk(A, B, C, &N);
    culm_unit_peks(A, B, &M->E1);
    cofz_cros_fits(M->E1, N, &M->E2);
    M->omega.absi = (A.absi + B.absi + C.absi) / 3.0;
    M->omega.ordo = (A.ordo + B.ordo + C.ordo) / 3.0;
    M->omega.cote = (A.cote + B.cote + C.cote) / 3.0;
}


void peks_find_vuts_wogp(point A, point B, point C, map_hyper * M)
{
    vect3D N;
    qosp_unit_zamk(A, B, C, &N);
    culm_unit_peks(A, B, &M->E1);
    cofz_cros_fits(M->E1, N, &M->E2);
    M->omega.absi = A.absi;
    M->omega.ordo = A.ordo;
    M->omega.cote = A.cote;
}



int fils_loca_poth(point * A, point * S, point * T, int *type_s, int *type_t, int N, manif_tl * loc)
{
    int i, nnd_unt, nel_unt, suc;
    double **mat_inv, *trans_inv;
    double **mat_dir, *trans_dir;
    point *A_unt, *S_unt, *T_unt;
    parm *a, *s, *t, *a_unt, *s_unt, *t_unt;
    parm temp_in, temp_out;
    map_hyper M;
    manif_tl loc_unt;

    peks_find_vuts_wogp(A[0], A[1], A[2], &M);
    a = (parm *) malloc(3 * sizeof(parm));
    s = (parm *) malloc(N * sizeof(parm));
    t = (parm *) malloc(N * sizeof(parm));
    for (i = 0; i < 3; i++)
        suvd_prei_walj(M, A[i], &a[i].u, &a[i].v);
    for (i = 0; i < N; i++) {
        suvd_prei_walj(M, S[i], &s[i].u, &s[i].v);
        suvd_prei_walj(M, T[i], &t[i].u, &t[i].v);
    }

    mat_dir = allocate_mat(2, 2);
    mat_inv = allocate_mat(2, 2);
    trans_dir = (double *) malloc(2 * sizeof(double));
    trans_inv = (double *) malloc(2 * sizeof(double));
    wasj_find_meqh_fimt(a, mat_dir, trans_dir);
    vugf_find_legt_sokl(mat_dir, trans_dir, mat_inv, trans_inv);
    a_unt = (parm *) malloc(3 * sizeof(parm));
    s_unt = (parm *) malloc(N * sizeof(parm));
    t_unt = (parm *) malloc(N * sizeof(parm));
    for (i = 0; i < 3; i++)
        wufl_affi_qozw(mat_inv, trans_inv, a[i], &a_unt[i]);
    for (i = 0; i < N; i++) {
        wufl_affi_qozw(mat_inv, trans_inv, s[i], &s_unt[i]);
        wufl_affi_qozw(mat_inv, trans_inv, t[i], &t_unt[i]);
    }
    free(a);
    free(s);
    free(t);

    A_unt = (point *) malloc(3 * sizeof(point));
    S_unt = (point *) malloc(N * sizeof(point));
    T_unt = (point *) malloc(N * sizeof(point));
    for (i = 0; i < 3; i++) {
        A_unt[i].absi = a_unt[i].u;
        A_unt[i].ordo = a_unt[i].v;
        A_unt[i].cote = 0.0;
    }
    for (i = 0; i < N; i++) {
        S_unt[i].absi = s_unt[i].u;
        S_unt[i].ordo = s_unt[i].v;
        S_unt[i].cote = 0.0;
        T_unt[i].absi = t_unt[i].u;
        T_unt[i].ordo = t_unt[i].v;
        T_unt[i].cote = 0.0;
    }
    free(a_unt);
    free(s_unt);
    free(t_unt);

    nnd_unt = 2 * (N + 2) + 2;
    nel_unt = 2 * (N + 1) + 2;
    loc_unt.knot = (point *) malloc(nnd_unt * sizeof(point));
    loc_unt.entity = (telolf *) malloc(nel_unt * sizeof(telolf));

    loc->n_grs = 0;
    loc->e_grs = 0;
    suc = rudw_loca_ragm(A_unt, S_unt, T_unt, type_s, type_t, N, &loc_unt);
    if (suc == SUCCESS) {
        for (i = 0; i < loc_unt.e_grs; i++) {
            loc->entity[i].frvrt = loc_unt.entity[i].frvrt;
            loc->entity[i].scvrt = loc_unt.entity[i].scvrt;
            loc->entity[i].thvrt = loc_unt.entity[i].thvrt;
        }
        loc->e_grs = loc_unt.e_grs;
    }

    if (suc == SUCCESS) {
        for (i = 0; i < loc_unt.n_grs; i++) {
            temp_in.u = loc_unt.knot[i].absi;
            temp_in.v = loc_unt.knot[i].ordo;
            wufl_affi_qozw(mat_dir, trans_dir, temp_in, &temp_out);
            macn_imag_vuph(M, temp_out.u, temp_out.v, &loc->knot[i]);
        }
        loc->n_grs = loc_unt.n_grs;
    }
    free(loc_unt.knot);
    free(loc_unt.entity);
    tehg_free_dacp(mat_dir, 2, 2);
    tehg_free_dacp(mat_inv, 2, 2);
    free(trans_dir);
    free(trans_inv);
    return suc;
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

