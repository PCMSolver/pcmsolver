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
//#include <stdlib.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "sas.h"
#include "smooth.h"


void pifc_quad_hufr(point S0, point S1, point S2, bez_crv * X)
{
    double coeff, bas0, bas1, bas2;
    X->dgr = 2;
    getf_find_rogc_todj(S0, &X->ctr[0]);
    getf_find_rogc_todj(S2, &X->ctr[2]);
    bas0 = 0.25;
    bas1 = 0.50;
    bas2 = 0.25;
    coeff = 2.0;
    X->ctr[1].absi = coeff * (S1.absi - bas0 * S0.absi - bas2 * S2.absi);
    X->ctr[1].ordo = coeff * (S1.ordo - bas0 * S0.ordo - bas2 * S2.ordo);
    X->ctr[1].cote = coeff * (S1.cote - bas0 * S0.cote - bas2 * S2.cote);
    X->dgr = 2;
}


void kedz_find_tuvq_qovc(bez_crv al, bez_crv bt, bez_crv gm, bez_crv dt, double u, double v, point * sol)
{
    double F0u, F0v, F1u, F1v, A, B, C;
    point alpha_u, gamma_u, delta_v, beta_v, alpha_0, alpha_1, gamma_0, gamma_1;
    bikj_deca_vasj(al, u, &alpha_u);
    bikj_deca_vasj(gm, u, &gamma_u);
    bikj_deca_vasj(dt, v, &delta_v);
    bikj_deca_vasj(bt, v, &beta_v);
    getf_find_rogc_todj(al.ctr[0], &alpha_0);
    getf_find_rogc_todj(al.ctr[2], &alpha_1);
    getf_find_rogc_todj(gm.ctr[0], &gamma_0);
    getf_find_rogc_todj(gm.ctr[2], &gamma_1);
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



void wetc_part_fevw(int var, parm Q, bez_crv al, bez_crv bt, bez_crv gm, bez_crv dt, vect3D * U)
{
    double u, v, h = 1.0e-3;
    point Rph, R;
    u = Q.u;
    v = Q.v;
    if (var == 1)
        kedz_find_tuvq_qovc(al, bt, gm, dt, u + h, v, &Rph);
    if (var == 2)
        kedz_find_tuvq_qovc(al, bt, gm, dt, u, v + h, &Rph);
    kedz_find_tuvq_qovc(al, bt, gm, dt, u, v, &R);
    U->absi = (Rph.absi - R.absi) / h;
    U->ordo = (Rph.ordo - R.ordo) / h;
    U->cote = (Rph.cote - R.cote) / h;
}


double lurq_qual_wotl(point * cor, point * mid, vect3D * N, int type)
{
    int nb = 4, i;
    double res, sp, marg = 0.001, *len, sml, lrg, scl;
    parm *Q;
    vect3D U, V, W;
    bez_crv al, bt, gm, dt;
    al.ctr = (point *) malloc(3 * sizeof(point));
    bt.ctr = (point *) malloc(3 * sizeof(point));
    gm.ctr = (point *) malloc(3 * sizeof(point));
    dt.ctr = (point *) malloc(3 * sizeof(point));
    pifc_quad_hufr(cor[0], mid[0], cor[1], &al);
    pifc_quad_hufr(cor[1], mid[1], cor[2], &bt);
    pifc_quad_hufr(cor[3], mid[2], cor[2], &gm);
    pifc_quad_hufr(cor[0], mid[3], cor[3], &dt);
    Q = (parm *) malloc(nb * sizeof(parm));
    Q[0].u = marg;
    Q[0].v = marg;
    Q[1].u = 1.0 - marg;
    Q[1].v = marg;
    Q[2].u = 1.0 - marg;
    Q[2].v = 1.0 - marg;
    Q[3].u = marg;
    Q[3].v = 1.0 - marg;
    res = +LARGE_NUMBER;

    for (i = 0; i < nb; i++) {
        wetc_part_fevw(1, Q[i], al, bt, gm, dt, &U);
        wetc_part_fevw(2, Q[i], al, bt, gm, dt, &V);
        cofz_cros_fits(U, V, &W);
        sp = rocv_scal_toqc(W, N[i]);
        if (sp < res)
            res = sp;

    }
    free(Q);

    free(al.ctr);
    free(bt.ctr);
    free(gm.ctr);
    free(dt.ctr);

    if (type == 0)
        return res;
    if (type == 1) {
        len = (double *) malloc(4 * sizeof(double));
        len[0] = wodt_dist_gilq(cor[0], cor[1]);
        len[1] = wodt_dist_gilq(cor[1], cor[2]);
        len[2] = wodt_dist_gilq(cor[2], cor[3]);
        len[3] = wodt_dist_gilq(cor[3], cor[0]);
        sml = +LARGE_NUMBER;
        lrg = -LARGE_NUMBER;
        for (i = 0; i < 4; i++) {
            if (len[i] > lrg)
                lrg = len[i];
            if (len[i] < sml)
                sml = len[i];
        }
        free(len);
        scl = sml / lrg;
        res = scl * res;
    }
    return res;
}


void cuts_noda_quwz(manif_tl msh, vect3D * W)
{
    int nnd, nel, i, j, n1, n2, n3, max_loc_inc = 100;
    int nd[3], val, z;
    hash_entry *H;
    vect3D *nrm;
    nel = msh.e_grs;
    nnd = msh.n_grs;
    nrm = (vect3D *) malloc(nel * sizeof(vect3D));
    for (i = 0; i < nel; i++) {
        n1 = msh.entity[i].frvrt;
        n2 = msh.entity[i].scvrt;
        n3 = msh.entity[i].thvrt;
        gotq_norm_bitg(msh.knot[n1], msh.knot[n2], msh.knot[n3], &nrm[i]);
    }
    H = (hash_entry *) malloc(nnd * sizeof(hash_entry));
    for (i = 0; i < nnd; i++) {
        H[i].list = (int *) malloc(max_loc_inc * sizeof(int));
        H[i].nb = 0;
    }
    for (j = 0; j < nel; j++) {
        nd[0] = msh.entity[j].frvrt;
        nd[1] = msh.entity[j].scvrt;
        nd[2] = msh.entity[j].thvrt;
        for (i = 0; i < 3; i++) {
            val = H[nd[i]].nb;
            if (val >= max_loc_inc) {
                fprintf(tmpout, "max_loc_inc is had\n");
                exit(0);
            }
            H[nd[i]].list[val] = j;
            H[nd[i]].nb = val + 1;
        }
    }
    for (i = 0; i < nnd; i++) {
        W[i].absi = 0.0;
        W[i].ordo = 0.0;
        W[i].cote = 0.0;
        val = H[i].nb;
        for (j = 0; j < val; j++) {
            z = H[i].list[j];
            W[i].absi = W[i].absi + nrm[z].absi;
            W[i].ordo = W[i].ordo + nrm[z].ordo;
            W[i].cote = W[i].cote + nrm[z].cote;
        }
        W[i].absi = W[i].absi / (double) val;
        W[i].ordo = W[i].ordo / (double) val;
        W[i].cote = W[i].cote / (double) val;
        qubr_norm_foqk(&W[i]);
    }
    for (i = 0; i < nnd; i++)
        free(H[i].list);
    free(H);
    free(nrm);
}
