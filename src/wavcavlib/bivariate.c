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
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "sas.h"
#include "smooth.h"


void copm_para_kesl(double *u, int L, double *v, int M, double *delta_u, double *delta_v)
{
    int i, k;
    for (i = 0; i <= L - 1; i++)
        delta_u[i] = u[i + 1] - u[i];
    for (k = 0; k <= M - 1; k++)
        delta_v[k] = v[k + 1] - v[k];
    delta_u[L] = 0.0;
    delta_v[M] = 0.0;
}


void patm_para_medh(double *delta, int N, double *sec)
{
    int i;
    for (i = 1; i < N; i++)
        sec[i] = delta[i] + delta[i - 1];
}


void vubq_find_qejg(int L, int M, double *delta_u, double *delta_v, double *sec_u, double *sec_v, double *Q_1, double *Q_2, double *D_1, double *D_2)
{
    int i, k;
    double den;
    for (i = 0; i <= L; i++) {
        if (i == 0)
            den = delta_u[i];
        else if (i == L)
            den = delta_u[i - 1];
        else
            den = sec_u[i];
        Q_1[i] = delta_u[i] / den;
        if (i > 0)
            Q_2[i] = delta_u[i - 1] / den;
        else
            Q_2[i] = 0.0;
    }
    for (k = 0; k <= M; k++) {
        if (k == 0)
            den = delta_v[k];
        else if (k == M)
            den = sec_v[k - 1];
        else
            den = sec_v[k];
        D_1[k] = delta_v[k] / den;
        if (k > 0)
            D_2[k] = delta_v[k - 1] / den;
        else
            D_2[k] = 0.0;
    }
}



void kogm_find_rinc(int L, int M, double *delta_u, double *delta_v, double *sec_u, double *A, double *B, double *C, double *D)
{
    int i;
    double den;
    for (i = 1; i <= L; i++) {
        if (i == 1)
            den = sec_u[i];
        else if (i == L)
            den = sec_u[i - 1];
        else
            den = sec_u[i] + delta_u[i - 2];
        if (i == L)
            A[i] = delta_u[i - 1] / den;
        else
            A[i] = sec_u[i] / den;
        if (i == 1)
            C[i] = 0.0;
        else
            C[i] = delta_u[i - 2] / den;
        B[i] = delta_u[i] / den;
        if (i == 1)
            D[i] = delta_u[i - 1] / den;
        else
            D[i] = sec_u[i - 1] / den;
    }
}



void qovk_find_nazv(double *delta_v, double *sec_v, int M, double *V_1, double *V_2, double *W_1, double *W_2)
{
    int k;
    double *DEN;
    DEN = (double *) malloc((M + 1) * sizeof(double));
    for (k = 1; k <= M; k++) {
        if (k == 1)
            DEN[k] = sec_v[k];
        else if (k == M)
            DEN[k] = sec_v[k - 1];
        else
            DEN[k] = sec_v[k] + delta_v[k - 2];
        V_1[k] = delta_v[k] / DEN[k];
        if (k == 1)
            W_1[k] = delta_v[k - 1] / DEN[k];
        else
            W_1[k] = sec_v[k - 1] / DEN[k];
    }
    for (k = 0; k <= M - 1; k++) {
        if (k == M - 1)
            V_2[k] = delta_v[k] / DEN[k + 1];
        else
            V_2[k] = sec_v[k + 1] / DEN[k + 1];
        if (k == 0)
            W_2[k] = 0.0;
        else
            W_2[k] = delta_v[k - 1] / DEN[k + 1];
    }
    free(DEN);
}


void cosg_asse_kens(int L, double *A, double *B, double *C, double *D, double *Q_1, double *Q_2, double **mat_lf)
{
    int i, j;
    double F1, F2, F3;
    for (i = 0; i <= L; i++)
        for (j = 0; j <= L; j++)
            mat_lf[i][j] = 0.0;
    mat_lf[0][0] = A[1];
    mat_lf[0][1] = C[1];
    for (i = 1; i <= L - 1; i++) {
        F1 = Q_1[i] * B[i];
        F2 = Q_1[i] * D[i] + Q_2[i] * A[i + 1];
        F3 = Q_2[i] * C[i + 1];
        mat_lf[i][i - 1] = F1;
        mat_lf[i][i] = F2;
        mat_lf[i][i + 1] = F3;
    }
    mat_lf[L][L - 1] = B[L];
    mat_lf[L][L] = D[L];
}



void fegn_asse_sevm(int M, double *V_1, double *V_2, double *W_1, double *W_2, double *D_1, double *D_2, double **mat_rg)
{
    int i, j, k;
    double G1, G2, G3;
    for (i = 0; i <= M; i++)
        for (j = 0; j <= M; j++)
            mat_rg[i][j] = 0.0;
    mat_rg[0][0] = V_2[0];
    mat_rg[1][0] = W_2[0];
    for (k = 1; k <= M - 1; k++) {
        G1 = D_1[k] * V_1[k];
        G2 = D_1[k] * W_1[k] + D_2[k] * V_2[k];
        G3 = D_2[k] * W_2[k];
        mat_rg[k - 1][k] = G1;
        mat_rg[k][k] = G2;
        mat_rg[k + 1][k] = G3;
    }
    mat_rg[M - 1][M] = V_1[M];
    mat_rg[M][M] = W_1[M];
}


void duhw_righ_vibq(double *u, double *v, int L, int M, double **mat_lf, double **mat_rg)
{
    double *delta_u, *delta_v;
    double *V_1, *V_2, *W_1, *W_2, *D_1, *D_2;
    double *A, *B, *C, *D, *Q_1, *Q_2, *sec_u, *sec_v;
    Q_1 = (double *) malloc((L + 1) * sizeof(double));
    Q_2 = (double *) malloc((L + 1) * sizeof(double));
    D_1 = (double *) malloc((M + 1) * sizeof(double));
    D_2 = (double *) malloc((M + 1) * sizeof(double));
    delta_u = (double *) malloc((L + 1) * sizeof(double));
    delta_v = (double *) malloc((M + 1) * sizeof(double));
    copm_para_kesl(u, L, v, M, delta_u, delta_v);
    sec_u = (double *) malloc((L + 1) * sizeof(double));
    sec_v = (double *) malloc((L + 1) * sizeof(double));
    patm_para_medh(delta_u, L, sec_u);
    patm_para_medh(delta_v, M, sec_v);
    vubq_find_qejg(L, M, delta_u, delta_v, sec_u, sec_v, Q_1, Q_2, D_1, D_2);
    V_1 = (double *) malloc((M + 1) * sizeof(double));
    V_2 = (double *) malloc(M * sizeof(double));
    W_1 = (double *) malloc((M + 1) * sizeof(double));
    W_2 = (double *) malloc(M * sizeof(double));
    qovk_find_nazv(delta_v, sec_v, M, V_1, V_2, W_1, W_2);
    A = (double *) malloc((L + 1) * sizeof(double));
    B = (double *) malloc((L + 1) * sizeof(double));
    C = (double *) malloc((L + 1) * sizeof(double));
    D = (double *) malloc((L + 1) * sizeof(double));
    kogm_find_rinc(L, M, delta_u, delta_v, sec_u, A, B, C, D);
    free(delta_u);
    free(delta_v);
    free(sec_u);
    free(sec_v);
    fegn_asse_sevm(M, V_1, V_2, W_1, W_2, D_1, D_2, mat_rg);
    free(V_1);
    free(V_2);
    free(W_1);
    free(W_2);
    cosg_asse_kens(L, A, B, C, D, Q_1, Q_2, mat_lf);
    free(A);
    free(B);
    free(C);
    free(D);
    free(Q_1);
    free(Q_2);
    free(D_1);
    free(D_2);
}


void hocv_corn_gohd(point omega, point A, point B, point C, point * X)
{
    double delta, lambda, coeff;
    vect3D S_1, S_2, W;
    delta = wodt_dist_gilq(omega, C);
    bofp_form_nukv(omega, A, &S_1);
    bofp_form_nukv(omega, B, &S_2);
    W.absi = S_1.absi + S_2.absi;
    W.ordo = S_1.ordo + S_2.ordo;
    W.cote = S_1.cote + S_2.cote;
    qubr_norm_foqk(&W);
    lambda = 1.0 / 3.0;
    coeff = lambda * delta;
    X->absi = omega.absi + coeff * W.absi;
    X->ordo = omega.ordo + coeff * W.ordo;
    X->cote = omega.cote + coeff * W.cote;
}


void hups_asse_runv(int L, int M, point ** P, point ** rhs)
{
    int i, j;
    double lambda, mu;
    point B_F, B_L;

    for (i = 1; i <= L - 1; i++)
        for (j = 1; j <= M - 1; j++)
            rhs[i][j] = P[i][j];

    lambda = 1.0 / 3.0;
    mu = 1.0 - lambda;
    for (j = 1; j <= M - 1; j++) {
        B_F.absi = mu * P[0][j].absi + lambda * P[1][j].absi;
        B_F.ordo = mu * P[0][j].ordo + lambda * P[1][j].ordo;
        B_F.cote = mu * P[0][j].cote + lambda * P[1][j].cote;

        B_L.absi = mu * P[L][j].absi + lambda * P[L - 1][j].absi;
        B_L.ordo = mu * P[L][j].ordo + lambda * P[L - 1][j].ordo;
        B_L.cote = mu * P[L][j].cote + lambda * P[L - 1][j].cote;

        getf_find_rogc_todj(B_F, &rhs[0][j]);
        getf_find_rogc_todj(B_L, &rhs[L][j]);
    }

    for (i = 1; i <= L - 1; i++) {
        B_F.absi = mu * P[i][0].absi + lambda * P[i][1].absi;
        B_F.ordo = mu * P[i][0].ordo + lambda * P[i][1].ordo;
        B_F.cote = mu * P[i][0].cote + lambda * P[i][1].cote;

        B_L.absi = mu * P[i][M].absi + lambda * P[i][M - 1].absi;
        B_L.ordo = mu * P[i][M].ordo + lambda * P[i][M - 1].ordo;
        B_L.cote = mu * P[i][M].cote + lambda * P[i][M - 1].cote;

        getf_find_rogc_todj(B_F, &rhs[i][0]);
        getf_find_rogc_todj(B_L, &rhs[i][M]);
    }

    hocv_corn_gohd(P[0][0], P[1][0], P[0][1], P[1][1], &rhs[0][0]);
    hocv_corn_gohd(P[0][M], P[0][M - 1], P[1][M], P[1][M - 1], &rhs[0][M]);
    hocv_corn_gohd(P[L][0], P[L][1], P[L - 1][0], P[L - 1][1], &rhs[L][0]);
    hocv_corn_gohd(P[L][M], P[L - 1][M], P[L][M - 1], P[L - 1][M - 1], &rhs[L][M]);
}


void nawm_solv_qerg(double **mat_lf, double **mat_rg, double **B, double **W, int L, int M)
{
    int i, j;
    double *s, *b, **S, **transp_Z;
    S = allocate_mat(L + 1, M + 1);

    s = (double *) malloc((L + 1) * sizeof(double));
    b = (double *) malloc((L + 1) * sizeof(double));
    for (j = 0; j <= M; j++) {
        for (i = 0; i <= L; i++)
            b[i] = B[i][j];
        poml_trid_fert(mat_lf, b, s, L + 1);
        for (i = 0; i <= L; i++)
            S[i][j] = s[i];
    }
    free(s);
    free(b);

    transp_Z = allocate_mat(M + 1, M + 1);
    for (i = 0; i <= M; i++)
        for (j = 0; j <= M; j++)
            transp_Z[i][j] = mat_rg[j][i];
    s = (double *) malloc((M + 1) * sizeof(double));
    b = (double *) malloc((M + 1) * sizeof(double));
    for (j = 0; j <= M; j++) {
        for (i = 0; i <= M; i++)
            b[i] = S[j][i];
        poml_trid_fert(transp_Z, b, s, M + 1);
        for (i = 0; i <= M; i++)
            W[j][i] = s[i];
    }

    free(s);
    free(b);
    tehg_free_dacp(transp_Z, M + 1, M + 1);
    tehg_free_dacp(S, L + 1, M + 1);
}


void ruhb_inte_vagr(int L, int M, double **mat_lf, double **mat_rg, point ** rhs, ns_surf * S)
{
    int i, j;
    double **B, **W;

    B = allocate_mat(L + 1, M + 1);
    W = allocate_mat(L + 1, M + 1);
    for (i = 0; i <= L; i++)
        for (j = 0; j <= M; j++)
            B[i][j] = rhs[i][j].absi;
    nawm_solv_qerg(mat_lf, mat_rg, B, W, L, M);
    for (i = 0; i <= L; i++)
        for (j = 0; j <= M; j++)
            S->d[i + 1][j + 1].absi = W[i][j];

    for (i = 0; i <= L; i++)
        for (j = 0; j <= M; j++)
            B[i][j] = rhs[i][j].ordo;
    nawm_solv_qerg(mat_lf, mat_rg, B, W, L, M);
    for (i = 0; i <= L; i++)
        for (j = 0; j <= M; j++)
            S->d[i + 1][j + 1].ordo = W[i][j];

    for (i = 0; i <= L; i++)
        for (j = 0; j <= M; j++)
            B[i][j] = rhs[i][j].cote;
    nawm_solv_qerg(mat_lf, mat_rg, B, W, L, M);
    for (i = 0; i <= L; i++)
        for (j = 0; j <= M; j++)
            S->d[i + 1][j + 1].cote = W[i][j];

    tehg_free_dacp(B, L + 1, M + 1);
    tehg_free_dacp(W, L + 1, M + 1);
}


void jepq_inte_lumr(int L, int M, ns_curv C_1, ns_curv C_2, ns_curv C_4, ns_curv C_3, ns_surf * S)
{
    int i, j;
    for (i = 0; i <= L + 2; i++) {
        getf_find_rogc_todj(C_1.d[i], &S->d[i][0]);
        getf_find_rogc_todj(C_2.d[i], &S->d[i][M + 2]);
    }
    for (j = 0; j <= M + 2; j++) {
        getf_find_rogc_todj(C_3.d[j], &S->d[0][j]);
        getf_find_rogc_todj(C_4.d[j], &S->d[L + 2][j]);
    }
}



void reck_bicu_bavh(double *u, double *v, point ** P, int L, int M, ns_curv C_1, ns_curv C_2, ns_curv C_3, ns_curv C_4, ns_surf * S)
{
    int i, k = 4;
    double **mat_lf, **mat_rg;
    point **rhs;

    mat_lf = allocate_mat(L + 1, L + 1);
    mat_rg = allocate_mat(M + 1, M + 1);
    rhs = allocate_mat_point(L + 1, M + 1);
    duhw_righ_vibq(u, v, L, M, mat_lf, mat_rg);
    hups_asse_runv(L, M, P, rhs);
    ruhb_inte_vagr(L, M, mat_lf, mat_rg, rhs, S);
    tehg_free_dacp(mat_lf, L + 1, L + 1);
    tehg_free_dacp(mat_rg, M + 1, M + 1);
    dirj_free_sukl(rhs, L + 1, M + 1);

    jepq_inte_lumr(L, M, C_1, C_2, C_4, C_3, S);
    S->nu = L + 2;
    S->nv = M + 2;
    S->ku = 4;
    S->kv = 4;

    for (i = 0; i < k; i++)
        S->frknt[i] = u[0];
    for (i = 1; i < L; i++)
        S->frknt[k + i - 1] = u[i];
    for (i = 1; i <= k; i++)
        S->frknt[L + 2 + i] = u[L];

    for (i = 0; i < k; i++)
        S->scknt[i] = u[0];
    for (i = 1; i < M; i++)
        S->scknt[k + i - 1] = v[i];
    for (i = 1; i <= k; i++)
        S->scknt[M + 2 + i] = v[M];

    S->prop1 = 0;
    S->prop2 = 0;
    S->prop3 = 1;
    S->prop4 = 0;
    S->prop5 = 0;
    S->u0 = u[0];
    S->u1 = u[L];
    S->v0 = v[0];
    S->v1 = v[M];
}


void qofm_bicu_bofl(double *u, double *v, point ** P, int L, int M, double **mat_lf, double **mat_rg, ns_curv C_1, ns_curv C_2, ns_curv C_3, ns_curv C_4, ns_surf * S)
{
    int i, k = 4;
    point **rhs;

    rhs = allocate_mat_point(L + 1, M + 1);
    hups_asse_runv(L, M, P, rhs);
    ruhb_inte_vagr(L, M, mat_lf, mat_rg, rhs, S);
    dirj_free_sukl(rhs, L + 1, M + 1);

    jepq_inte_lumr(L, M, C_1, C_2, C_4, C_3, S);
    S->nu = L + 2;
    S->nv = M + 2;
    S->ku = 4;
    S->kv = 4;

    for (i = 0; i < k; i++)
        S->frknt[i] = u[0];
    for (i = 1; i < L; i++)
        S->frknt[k + i - 1] = u[i];
    for (i = 1; i <= k; i++)
        S->frknt[L + 2 + i] = u[L];

    for (i = 0; i < k; i++)
        S->scknt[i] = u[0];
    for (i = 1; i < M; i++)
        S->scknt[k + i - 1] = v[i];
    for (i = 1; i <= k; i++)
        S->scknt[M + 2 + i] = v[M];

    S->prop1 = 0;
    S->prop2 = 0;
    S->prop3 = 1;
    S->prop4 = 0;
    S->prop5 = 0;
    S->u0 = u[0];
    S->u1 = u[L];
    S->v0 = v[0];
    S->v1 = v[M];
}
