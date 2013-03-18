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
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"



double vikr_dist_rugc(point A, point B, point C, point X)
{
    double res;
    vect3D AB, AC, W, AX;
    bofp_form_nukv(A, B, &AB);
    bofp_form_nukv(A, C, &AC);
    cofz_cros_fits(AB, AC, &W);
    qubr_norm_foqk(&W);
    bofp_form_nukv(A, X, &AX);
    res = rocv_scal_toqc(AX, W);
    return fabs(res);
}



void rehp_find_deng(point * P, int n, int id_0, int id_1, int id_2, int *idx)
{
    int i, i_sp, dummy, ts;
    double lg, ds;
    idx[0] = id_0;
    idx[1] = id_1;
    idx[2] = id_2;

    lg = 0.0;
    for (i = 0; i < n; i++) {
        ts = gonl_arra_govj(idx, 3, i, &dummy);
        if (ts == 0) {
            ds = vikr_dist_rugc(P[idx[0]], P[idx[1]], P[idx[2]], P[i]);
            if (ds > lg) {
                lg = ds;
                i_sp = i;
            }
        }
    }
    idx[3] = i_sp;
}



void qumg_circ_tezl(point a, point b, point c, point d, double *rad, point * center)
{
    double den, s_ab, s_ac, s_ad, nrm;
    vect3D ab, ac, ad, c1, c2, c3, num;
    bofp_form_nukv(a, b, &ab);
    bofp_form_nukv(a, c, &ac);
    bofp_form_nukv(a, d, &ad);
    cofz_cros_fits(ac, ad, &c1);
    cofz_cros_fits(ab, ac, &c2);
    cofz_cros_fits(ad, ab, &c3);
    den = 2.0 * rocv_scal_toqc(ab, c1);
    s_ab = rocv_scal_toqc(ab, ab);
    s_ac = rocv_scal_toqc(ac, ac);
    s_ad = rocv_scal_toqc(ad, ad);
    num.absi = s_ad * c2.absi + s_ac * c3.absi + s_ab * c1.absi;
    num.ordo = s_ad * c2.ordo + s_ac * c3.ordo + s_ab * c1.ordo;
    num.cote = s_ad * c2.cote + s_ac * c3.cote + s_ab * c1.cote;
    center->absi = a.absi + (num.absi / den);
    center->ordo = a.ordo + (num.ordo / den);
    center->cote = a.cote + (num.cote / den);
    nrm = biqh_norm_dapf(num);
    *rad = nrm / fabs(den);
}


int qojc_reco_fadc(point * P, int n, int id_0, int id_1, int id_2, double eps, point * center, double *rad)
{
    int *idx, i, s, dummy, res, ts, nb;
    double dist, rel_eps, RAD, total_err;
    double xmi, xma, ymi, yma, zmi, zma, h_x, h_y, h_z, sr[3], surf;
    point *Q;

    idx = (int *) malloc(4 * sizeof(int));
    rehp_find_deng(P, n, id_0, id_1, id_2, idx);
    Q = (point *) malloc(4 * sizeof(point));
    for (i = 0; i < 4; i++) {
        s = idx[i];
        getf_find_rogc_todj(P[s], &Q[i]);
    }
    qumg_circ_tezl(Q[0], Q[1], Q[2], Q[3], rad, center);

    RAD = *rad;
    rel_eps = RAD * eps;
    total_err = 0.0;
    nb = 0;
    for (i = 0; i < n; i++) {
        ts = gonl_arra_govj(idx, 4, i, &dummy);
        if (ts == 0) {
            dist = wodt_dist_gilq(*center, P[i]);
            total_err = total_err + fabs(dist - RAD);
            nb++;
        }
    }
    total_err = total_err / (double) nb;
    res = 1;
    if (total_err > rel_eps)
        res = 0;
    free(idx);
    free(Q);
    if (res == 1) {
        homs_boun_gosm(P, n, &xmi, &xma, &ymi, &yma, &zmi, &zma);
        h_x = xma - xmi;
        h_y = yma - ymi;
        h_z = zma - zmi;
        sr[0] = h_x * h_y;
        sr[1] = h_x * h_z;
        sr[2] = h_y * h_z;
        surf = 0.0;
        for (i = 0; i < 3; i++)
            if (sr[i] > surf)
                surf = sr[i];
        if (surf < 0.05 * RAD)
            res = 0;
    }
    return res;
}



void tefq_find_lebg_nehj(parm * A, parm X, baryc2D * q)
{
    int i;
    double **mat, *rhs, *lambda;
    mat = (double **) malloc(3 * sizeof(double *));
    for (i = 0; i < 3; i++)
        mat[i] = (double *) malloc(3 * sizeof(double));
    for (i = 0; i < 3; i++) {
        mat[0][i] = 1.0;
        mat[1][i] = A[i].u;
        mat[2][i] = A[i].v;
    }
    rhs = (double *) malloc(3 * sizeof(double));
    rhs[0] = 1.0;
    rhs[1] = X.u;
    rhs[2] = X.v;
    lambda = (double *) malloc(3 * sizeof(double));
    javs_find_gokm_pedt(mat, rhs, lambda, 3);
    q->lambda1 = lambda[0];
    q->lambda2 = lambda[1];
    q->lambda3 = lambda[2];
    free(lambda);
    free(rhs);
    for (i = 0; i < 3; i++)
        free(mat[i]);
    free(mat);
}



void demn_tria_wusg(manif_ro msh_tri, pt_cnv PC, manif_tl * msh)
{
    int nnd, nel, i, n1, n2, n3;
    double sp;
    point A, B, C, concavity, G;
    vect3D U, V;
    parm Q, *a;
    baryc2D q;

    a = (parm *) malloc(3 * sizeof(parm));
    a[0].u = 0.0;
    a[0].v = 0.0;
    a[1].u = 1.0;
    a[1].v = 0.0;
    a[2].u = 0.0;
    a[2].v = 1.0;
    nnd = msh_tri.n_grs;
    for (i = 0; i < nnd; i++) {
        if (PC.sitn == 0)
            rows_eval_qusg(PC, msh_tri.knot[i], &msh->knot[i]);
        if (PC.sitn == 2) {
            tefq_find_lebg_nehj(a, msh_tri.knot[i], &q);
            lekr_find_sujr_niqc(q, PC.cc, &Q);
            rows_eval_qusg(PC, Q, &msh->knot[i]);
        }
    }
    free(a);
    msh->n_grs = nnd;

    nel = msh_tri.e_grs;
    for (i = 0; i < nel; i++) {
        msh->entity[i].frvrt = msh_tri.entity[i].frvrt;
        msh->entity[i].scvrt = msh_tri.entity[i].scvrt;
        msh->entity[i].thvrt = msh_tri.entity[i].thvrt;
    }
    msh->e_grs = nel;

    getf_find_rogc_todj(PC.alpha.begn, &A);
    getf_find_rogc_todj(PC.beta.begn, &B);
    getf_find_rogc_todj(PC.gamma.begn, &C);
    concavity.absi = (A.absi + B.absi + C.absi) / 3.0;
    concavity.ordo = (A.ordo + B.ordo + C.ordo) / 3.0;
    concavity.cote = (A.cote + B.cote + C.cote) / 3.0;

    n1 = msh->entity[0].frvrt;
    n2 = msh->entity[0].scvrt;
    n3 = msh->entity[0].thvrt;
    G.absi = (msh->knot[n1].absi + msh->knot[n2].absi + msh->knot[n3].absi) / 3.0;
    G.ordo = (msh->knot[n1].ordo + msh->knot[n2].ordo + msh->knot[n3].ordo) / 3.0;
    G.cote = (msh->knot[n1].cote + msh->knot[n2].cote + msh->knot[n3].cote) / 3.0;
    culm_unit_peks(G, concavity, &U);
    gotq_norm_bitg(msh->knot[n1], msh->knot[n2], msh->knot[n3], &V);
    sp = rocv_scal_toqc(U, V);
    if (sp < 0.0) {
        for (i = 0; i < nel; i++) {
            n1 = msh->entity[i].frvrt;
            n2 = msh->entity[i].scvrt;
            n3 = msh->entity[i].thvrt;
            msh->entity[i].frvrt = n1;
            msh->entity[i].scvrt = n3;
            msh->entity[i].thvrt = n2;
        }
    }
}



void kotr_find_tujg(manif_ro msh_tri, int id_0, int id_1, int id_2, pt_cnv PC, conc_model * CM)
{
    int i, nel, nnd, n1, n2, n3, ts;
    double eps = 0.1, rad, sp;
    point center, *A, CT;
    vect3D U1, U2;
    demn_tria_wusg(msh_tri, PC, &CM->msh);
    nnd = CM->msh.n_grs;
    nel = CM->msh.e_grs;
    homs_boun_gosm(CM->msh.knot, nnd, &CM->BB.x_min, &CM->BB.x_max, &CM->BB.y_min, &CM->BB.y_max, &CM->BB.z_min, &CM->BB.z_max);

    for (i = 0; i < nel; i++) {
        n1 = CM->msh.entity[i].frvrt;
        n2 = CM->msh.entity[i].scvrt;
        n3 = CM->msh.entity[i].thvrt;
        CM->G[i].absi = (CM->msh.knot[n1].absi + CM->msh.knot[n2].absi + CM->msh.knot[n3].absi) / 3.0;
        CM->G[i].ordo = (CM->msh.knot[n1].ordo + CM->msh.knot[n2].ordo + CM->msh.knot[n3].ordo) / 3.0;
        CM->G[i].cote = (CM->msh.knot[n1].cote + CM->msh.knot[n2].cote + CM->msh.knot[n3].cote) / 3.0;
        gotq_norm_bitg(CM->msh.knot[n1], CM->msh.knot[n2], CM->msh.knot[n3], &CM->nrm[i]);
    }

    ts = qojc_reco_fadc(CM->msh.knot, nnd, id_0, id_1, id_2, eps, &center, &rad);
    if (ts == 0) {
        fprintf(tmpout, "Unable to find supporting sphere\n");
        exit(0);
    }
    getf_find_rogc_todj(center, &CM->S.zent);
    CM->S.rad = rad;

    A = (point *) malloc(3 * sizeof(point));
    getf_find_rogc_todj(PC.alpha.begn, &A[0]);
    getf_find_rogc_todj(PC.beta.begn, &A[1]);
    getf_find_rogc_todj(PC.gamma.begn, &A[2]);
    CT.absi = (A[0].absi + A[1].absi + CM->S.zent.absi) / 3.0;
    CT.ordo = (A[0].ordo + A[1].ordo + CM->S.zent.ordo) / 3.0;
    CT.cote = (A[0].cote + A[1].cote + CM->S.zent.cote) / 3.0;
    gotq_norm_bitg(CM->S.zent, A[0], A[1], &U1);
    bofp_form_nukv(CT, A[2], &U2);
    sp = rocv_scal_toqc(U1, U2);
    if (sp >= 0) {
        getf_find_rogc_todj(A[0], &CM->cr0);
        getf_find_rogc_todj(A[1], &CM->cr1);
        getf_find_rogc_todj(A[2], &CM->cr2);
    } else {
        getf_find_rogc_todj(A[0], &CM->cr0);
        getf_find_rogc_todj(A[2], &CM->cr1);
        getf_find_rogc_todj(A[1], &CM->cr2);
    }
    free(A);
}



int sarq_test_durp(point X, manif_tl msh, point * G, vect3D * nrm)
{
    int nel, i, res;
    double sp;
    vect3D U;
    res = 1;
    nel = msh.e_grs;
    for (i = 0; i < nel; i++) {
        culm_unit_peks(G[i], X, &U);
        sp = rocv_scal_toqc(U, nrm[i]);
        if (sp <= 0.0) {
            res = 0;
            break;
        }
    }
    return res;
}


int zejn_insi_jurg(point X, conc_model CM)
{
    int res = 1;
    double sp;
    point CT;
    vect3D U1, U2;
    CT.absi = (CM.cr0.absi + CM.cr1.absi + CM.S.zent.absi) / 3.0;
    CT.ordo = (CM.cr0.ordo + CM.cr1.ordo + CM.S.zent.ordo) / 3.0;
    CT.cote = (CM.cr0.cote + CM.cr1.cote + CM.S.zent.cote) / 3.0;
    gotq_norm_bitg(CM.S.zent, CM.cr0, CM.cr1, &U1);
    bofp_form_nukv(CT, X, &U2);
    sp = rocv_scal_toqc(U1, U2);
    if (sp <= 0)
        res = 0;

    if (res == 1) {
        CT.absi = (CM.cr1.absi + CM.cr2.absi + CM.S.zent.absi) / 3.0;
        CT.ordo = (CM.cr1.ordo + CM.cr2.ordo + CM.S.zent.ordo) / 3.0;
        CT.cote = (CM.cr1.cote + CM.cr2.cote + CM.S.zent.cote) / 3.0;
        gotq_norm_bitg(CM.S.zent, CM.cr1, CM.cr2, &U1);
        bofp_form_nukv(CT, X, &U2);
        sp = rocv_scal_toqc(U1, U2);
        if (sp <= 0)
            res = 0;
    }

    if (res == 1) {
        CT.absi = (CM.cr2.absi + CM.cr0.absi + CM.S.zent.absi) / 3.0;
        CT.ordo = (CM.cr2.ordo + CM.cr0.ordo + CM.S.zent.ordo) / 3.0;
        CT.cote = (CM.cr2.cote + CM.cr0.cote + CM.S.zent.cote) / 3.0;
        gotq_norm_bitg(CM.S.zent, CM.cr2, CM.cr0, &U1);
        bofp_form_nukv(CT, X, &U2);
        sp = rocv_scal_toqc(U1, U2);
        if (sp <= 0)
            res = 0;
    }
    return res;
}


int sunk_inte_qign(conc_model CM1, conc_model CM2)
{
    int nnd1, nnd2, i, ts, res, tr;
    nnd2 = CM2.msh.n_grs;
    res = 0;
    for (i = 0; i < nnd2; i++) {
        ts = sarq_test_durp(CM2.msh.knot[i], CM1.msh, CM1.G, CM1.nrm);
        if (ts == 1) {
            tr = zejn_insi_jurg(CM2.msh.knot[i], CM1);
            if (tr == 1) {
                res = 1;
                break;
            }
        }
    }

    if (res == 0) {
        nnd1 = CM1.msh.n_grs;
        for (i = 0; i < nnd1; i++) {
            ts = sarq_test_durp(CM1.msh.knot[i], CM2.msh, CM2.G, CM2.nrm);
            if (ts == 1) {
                tr = zejn_insi_jurg(CM1.msh.knot[i], CM2);
                if (tr == 1) {
                    res = 1;
                    break;
                }
            }
        }
    }
    return res;
}


int ticj_inte_ludj(conc_model CM1, conc_model CM2)
{
    int ts_interf, res;
    double dis, D, rad1, rad2, marg = 0.1;
    point cent1, cent2;
    ts_interf = lafc_boun_gusd(CM1.BB, CM2.BB);
    if (ts_interf == 0)
        return 0;
    getf_find_rogc_todj(CM1.S.zent, &cent1);
    getf_find_rogc_todj(CM2.S.zent, &cent2);
    dis = wodt_dist_gilq(cent1, cent2);
    rad1 = CM1.S.rad;
    rad2 = CM2.S.rad;
    D = rad1 + rad2;
    if (dis > D + marg)
        return 0;

    res = sunk_inte_qign(CM1, CM2);
    return res;
}


void ritc_find_cugz_jotw(manif_ro * msh)
{
    int nel, i, n1, n2, n3;
    vect2D U, V;
    vect3D U3D, V3D, Z;
    nel = msh->e_grs;
    for (i = 0; i < nel; i++) {
        n1 = msh->entity[i].frvrt;
        n2 = msh->entity[i].scvrt;
        n3 = msh->entity[i].thvrt;
        cuwl_unit_pist(msh->knot[n1], msh->knot[n2], &U);
        cuwl_unit_pist(msh->knot[n1], msh->knot[n3], &V);
        U3D.absi = U.u;
        V3D.absi = V.u;
        U3D.ordo = U.v;
        V3D.ordo = V.v;
        U3D.cote = 0.0;
        V3D.cote = 0.0;
        cofz_cros_fits(U3D, V3D, &Z);
        if (Z.cote < 0.0) {
            msh->entity[i].scvrt = n3;
            msh->entity[i].thvrt = n2;
        }
    }
}



void fiwh_simm_tucs(int N, manif_ro * msh, int *id_0, int *id_1, int *id_2)
{
    int nnd, nel, k, i, ts, id, *map;
    int n1, n2, n3, cr_a, cr_b, cr_c;
    double eps = 1.0e-9;
    parm A, B, C;
    manif_ro temp;
    A.u = 0.0;
    A.v = 0.0;
    B.u = 1.0;
    B.v = 0.0;
    C.u = 0.0;
    C.v = 1.0;
    nnd = 3 * N * N;
    nel = 6 * (N - 1) * (N - 1);
    temp.knot = (parm *) malloc(nnd * sizeof(parm));
    temp.entity = (telolf *) malloc(nel * sizeof(telolf));
    guqw_simp_howc(N, A, B, C, &temp, &cr_a, &cr_b, &cr_c);

    map = (int *) malloc(nnd * sizeof(int));
    k = 0;
    for (i = 0; i < nnd; i++) {
        ts = keld_arra_kefg(msh->knot, k, temp.knot[i], eps, &id);
        if (ts == 0) {
            cunl_find_qedf_rewn(temp.knot[i], &msh->knot[k]);
            map[i] = k;
            k++;
        } else
            map[i] = id;
        if (i == cr_a)
            *id_0 = map[i];
        if (i == cr_b)
            *id_1 = map[i];
        if (i == cr_c)
            *id_2 = map[i];
    }
    msh->n_grs = k;

    for (i = 0; i < nel; i++) {
        n1 = temp.entity[i].frvrt;
        n2 = temp.entity[i].scvrt;
        n3 = temp.entity[i].thvrt;
        msh->entity[i].frvrt = map[n1];
        msh->entity[i].scvrt = map[n2];
        msh->entity[i].thvrt = map[n3];
    }
    msh->e_grs = nel;
    ritc_find_cugz_jotw(msh);
    free(temp.entity);
    free(temp.knot);
    free(map);
}
