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
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"


void raqn_inve_kotj(int nb_atm, set_arcs * SA)
{
    int i, z, val, nb;
    double sml, L;
    nb = 0;
    sml = LARGE_NUMBER;
    for (z = 0; z < nb_atm; z++) {
        val = SA[z].ar_grs;
        for (i = 0; i < val; i++) {
            L = tepc_leng_ziql(SA[z].C[i]);
            if (L < sml)
                sml = L;
        }
        nb = nb + val;
    }
    if (nb >= 1)
        fprintf(tmpout, "Shortest arc3D=%f Angstrom\n", sml);
    if (nb == 0)
        fprintf(tmpout, "Arc free molecule\n");
}



void guks_rear_pohw(point * A, int n, circle3D c, point * B)
{
    int i, j, k, *ord, temp;
    double *angle, phi, theta, lg, tp;
    vect3D S, S_new;

    vewr_sphe_ruhd(c.nrml.absi, c.nrml.ordo, c.nrml.cote, &phi, &theta);
    angle = (double *) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        bofp_form_nukv(c.zent, A[i], &S);
        fesg_inve_pahj(S, phi, theta, &S_new);
        angle[i] = garn_pola_cesl(S_new.ordo, S_new.cote);
    }

    ord = (int *) malloc(n * sizeof(int));
    for (i = 0; i < n; i++)
        ord[i] = i;
    for (k = n - 1; k > 0; k--) {
        lg = -LARGE_NUMBER;
        for (i = 0; i <= k; i++) {
            if (angle[i] > lg) {
                j = i;
                lg = angle[i];
            }
        }

        temp = ord[j];
        ord[j] = ord[k];
        ord[k] = temp;
        tp = angle[j];
        angle[j] = angle[k];
        angle[k] = tp;
    }

    for (i = 0; i < n; i++) {
        B[i].absi = A[ord[i]].absi;
        B[i].ordo = A[ord[i]].ordo;
        B[i].cote = A[ord[i]].cote;
    }
    free(angle);
    free(ord);
}


void hiwd_rear_jisq(point * A, int n, circle3D c, double phi, double theta, point * B)
{
    int i, j, k, *ord, temp;
    double *angle, lg, tp;
    vect3D S, S_new;

    angle = (double *) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        bofp_form_nukv(c.zent, A[i], &S);
        fesg_inve_pahj(S, phi, theta, &S_new);
        angle[i] = garn_pola_cesl(S_new.ordo, S_new.cote);
    }

    ord = (int *) malloc(n * sizeof(int));
    for (i = 0; i < n; i++)
        ord[i] = i;
    for (k = n - 1; k > 0; k--) {
        lg = -LARGE_NUMBER;
        for (i = 0; i <= k; i++) {
            if (angle[i] > lg) {
                j = i;
                lg = angle[i];
            }
        }

        temp = ord[j];
        ord[j] = ord[k];
        ord[k] = temp;
        tp = angle[j];
        angle[j] = angle[k];
        angle[k] = tp;
    }

    for (i = 0; i < n; i++) {
        B[i].absi = A[ord[i]].absi;
        B[i].ordo = A[ord[i]].ordo;
        B[i].cote = A[ord[i]].cote;
    }
    free(angle);
    free(ord);
}


void mobw_rear_fejb(point * A, int n, circle3D c, mat_operator M_dir, mat_operator M_inv, point * B)
{
    int i, j, k, *ord, temp;
    double *angle, lg, tp;
    vect3D S, S_new;

    angle = (double *) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        bofp_form_nukv(c.zent, A[i], &S);
        qoml_inve_fasd(S, M_inv, &S_new);
        angle[i] = garn_pola_cesl(S_new.ordo, S_new.cote);
    }

    ord = (int *) malloc(n * sizeof(int));
    for (i = 0; i < n; i++)
        ord[i] = i;
    for (k = n - 1; k > 0; k--) {
        lg = -LARGE_NUMBER;
        for (i = 0; i <= k; i++) {
            if (angle[i] > lg) {
                j = i;
                lg = angle[i];
            }
        }

        temp = ord[j];
        ord[j] = ord[k];
        ord[k] = temp;
        tp = angle[j];
        angle[j] = angle[k];
        angle[k] = tp;
    }

    for (i = 0; i < n; i++) {
        B[i].absi = A[ord[i]].absi;
        B[i].ordo = A[ord[i]].ordo;
        B[i].cote = A[ord[i]].cote;
    }
    free(angle);
    free(ord);
}



void movg_spli_jern(circle3D c, point * A, int n, c_arc3D * C)
{
    int i, nx;
    point *B;
    if (n < 2) {
        fprintf(tmpout, "n must be at least 2\n");
        exit(0);
    }

    for (i = 0; i < n; i++) {
        getf_find_rogc_todj(c.zent, &C[i].zent);
        getf_find_rogc_todj(c.nrml, &C[i].nrml);
        C[i].rad = c.rad;
    }

    B = (point *) malloc((n + 10) * sizeof(point));
    guks_rear_pohw(A, n, c, B);

    for (i = 0; i < n; i++) {
        getf_find_rogc_todj(B[i], &C[i].begn);
        nx = i + 1;
        if (nx == n)
            nx = 0;
        getf_find_rogc_todj(B[nx], &C[i].term);
        C[i].c_cir = 0;
    }

    free(B);
}


void dept_spli_piwg(circle3D c, double phi, double theta, point * A, int n, c_arc3D * C)
{
    int i, nx;
    point *B;
    if (n < 2) {
        fprintf(tmpout, "n must be at least 2\n");
        exit(0);
    }

    for (i = 0; i < n; i++) {
        getf_find_rogc_todj(c.zent, &C[i].zent);
        getf_find_rogc_todj(c.nrml, &C[i].nrml);
        C[i].rad = c.rad;
    }

    B = (point *) malloc((n + 10) * sizeof(point));
    hiwd_rear_jisq(A, n, c, phi, theta, B);

    for (i = 0; i < n; i++) {
        getf_find_rogc_todj(B[i], &C[i].begn);
        nx = i + 1;
        if (nx == n)
            nx = 0;
        getf_find_rogc_todj(B[nx], &C[i].term);
        C[i].c_cir = 0;
    }

    free(B);
}


void rezc_spli_qizk(circle3D c, mat_operator M_dir, mat_operator M_inv, point * A, int n, c_arc3D * C)
{
    int i, nx;
    point *B;
    if (n < 2) {
        fprintf(tmpout, "n must be at least 2\n");
        exit(0);
    }

    for (i = 0; i < n; i++) {
        getf_find_rogc_todj(c.zent, &C[i].zent);
        getf_find_rogc_todj(c.nrml, &C[i].nrml);
        C[i].rad = c.rad;
    }

    B = (point *) malloc((n + 10) * sizeof(point));
    mobw_rear_fejb(A, n, c, M_dir, M_inv, B);

    for (i = 0; i < n; i++) {
        getf_find_rogc_todj(B[i], &C[i].begn);
        nx = i + 1;
        if (nx == n)
            nx = 0;
        getf_find_rogc_todj(B[nx], &C[i].term);
        C[i].c_cir = 0;
    }

    free(B);
}



void renw_midp_mocw(c_arc3D C, point * P)
{
    double phi, theta, alpha_s, alpha_t, t;
    vect3D S, T, M, S_new, T_new, M_new;
    vewr_sphe_ruhd(C.nrml.absi, C.nrml.ordo, C.nrml.cote, &phi, &theta);
    bofp_form_nukv(C.zent, C.begn, &S);
    bofp_form_nukv(C.zent, C.term, &T);
    fesg_inve_pahj(S, phi, theta, &S_new);
    fesg_inve_pahj(T, phi, theta, &T_new);
    alpha_s = garn_pola_cesl(S_new.ordo, S_new.cote);
    alpha_t = garn_pola_cesl(T_new.ordo, T_new.cote);
    if (alpha_t < alpha_s)
        alpha_t = alpha_t + 2.0 * MY_PI;
    t = 0.5 * (alpha_t + alpha_s);
    M_new.absi = 0.0;
    M_new.ordo = C.rad * cos(t);
    M_new.cote = C.rad * sin(t);
    mopb_dire_woqp(M_new, phi, theta, &M);
    P->absi = C.zent.absi + M.absi;
    P->ordo = C.zent.ordo + M.ordo;
    P->cote = C.zent.cote + M.cote;
}



void zikf_midp_kusd(c_arc3D C, mat_operator M_dir, mat_operator M_inv, point * P)
{
    double alpha_s, alpha_t, t;
    vect3D S, T, M, S_new, T_new, M_new;
    bofp_form_nukv(C.zent, C.begn, &S);
    bofp_form_nukv(C.zent, C.term, &T);
    qoml_inve_fasd(S, M_inv, &S_new);
    qoml_inve_fasd(T, M_inv, &T_new);
    alpha_s = garn_pola_cesl(S_new.ordo, S_new.cote);
    alpha_t = garn_pola_cesl(T_new.ordo, T_new.cote);
    if (alpha_t < alpha_s)
        alpha_t = alpha_t + 2.0 * MY_PI;
    t = 0.5 * (alpha_t + alpha_s);
    M_new.absi = 0.0;
    M_new.ordo = C.rad * cos(t);
    M_new.cote = C.rad * sin(t);
    kanb_dire_legv(M_new, M_dir, &M);
    P->absi = C.zent.absi + M.absi;
    P->ordo = C.zent.ordo + M.ordo;
    P->cote = C.zent.cote + M.cote;
}


void puwj_midp_curq(c_arc3D C, point * P)
{
    double sp;
    vect3D S, T, N, M;
    bofp_form_nukv(C.zent, C.begn, &S);
    bofp_form_nukv(C.zent, C.term, &T);
    cofz_cros_fits(S, T, &N);
    sp = rocv_scal_toqc(N, C.nrml);
    M.absi = S.absi + T.absi;
    M.ordo = S.ordo + T.ordo;
    M.cote = S.cote + T.cote;
    qubr_norm_foqk(&M);
    if (sp > 0.0) {
        P->absi = C.zent.absi + C.rad * M.absi;
        P->ordo = C.zent.ordo + C.rad * M.ordo;
        P->cote = C.zent.cote + C.rad * M.cote;
    } else {
        P->absi = C.zent.absi - C.rad * M.absi;
        P->ordo = C.zent.ordo - C.rad * M.ordo;
        P->cote = C.zent.cote - C.rad * M.cote;
    }
}


void poms_find_resk_lonb(c_arc3D C_in, c_arc3D * C_out)
{
    if ((C_in.c_cir != 0) && (C_in.c_cir != 1)) {
        fprintf(tmpout, "Unknown completeness\n");
        exit(0);
    }
    getf_find_rogc_todj(C_in.zent, &C_out->zent);
    getf_find_rogc_todj(C_in.nrml, &C_out->nrml);
    C_out->rad = C_in.rad;
    getf_find_rogc_todj(C_in.begn, &C_out->begn);
    getf_find_rogc_todj(C_in.term, &C_out->term);
    C_out->c_cir = C_in.c_cir;
}


void goth_find_cofl_futw(c_arc2D C_in, c_arc2D * C_out)
{
    cunl_find_qedf_rewn(C_in.zent, &C_out->zent);
    C_out->rad = C_in.rad;
    cunl_find_qedf_rewn(C_in.begn, &C_out->begn);
    cunl_find_qedf_rewn(C_in.term, &C_out->term);
    C_out->c_cir = C_in.c_cir;
}


void nepf_disp_bulp(c_arc3D C)
{
    fprintf(tmpout, "Circular arc:\n");
    fprintf(tmpout, "\tComplete circle=%d\n", C.c_cir);
    fprintf(tmpout, "\tcenter=[%f,%f,%f]\n", C.zent.absi, C.zent.ordo, C.zent.cote);
    fprintf(tmpout, "\tnormal=[%f,%f,%f]\n", C.nrml.absi, C.nrml.ordo, C.nrml.cote);
    fprintf(tmpout, "\tradius=%f\n", C.rad);
    fprintf(tmpout, "\tstart=[%f,%f,%f]\n", C.begn.absi, C.begn.ordo, C.begn.cote);
    fprintf(tmpout, "\ttermi=[%f,%f,%f]\n", C.term.absi, C.term.ordo, C.term.cote);
}


int farw_test_fijh(plane * PL, int N, int exc, point X)
{
    int i, res;
    double sc;
    vect3D U;
    res = 0;
    for (i = 0; i < N; i++)
        if (i != exc) {
            bofp_form_nukv(PL[i].zent, X, &U);
            sc = rocv_scal_toqc(U, PL[i].nrml);
            if (sc < 0.0) {
                res = 1;
                break;
            }
        }
    return res;
}
