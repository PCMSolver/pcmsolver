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
#include "splinemol.h"
#include "pln_sph.h"
#include "coarsequad.h"



void kowv_loca_colb(vect3D N, vect3D * U, vect3D * V)
{
    point t;
    double a, b, c, w;
    a = N.absi;
    b = N.ordo;
    c = N.cote;

    w = a * a + b * b;
    if (w < 0.00001) {
        t.absi = +1.0;
        t.ordo = 0.0;
        t.cote = 0.0;
    } else if (fabs(c) < 0.00001) {
        t.absi = -b;
        t.ordo = a;
        t.cote = 0.0;
    } else {
        t.absi = a;
        t.ordo = b;
        t.cote = -(a * a + b * b) / c;
    }

    U->absi = t.absi;
    U->ordo = t.ordo;
    U->cote = t.cote;
    qubr_norm_foqk(U);

    cofz_cros_fits(N, *U, V);
    qubr_norm_foqk(V);
}



void gufw_plan_mesv(vect3D U, vect3D V, point X, parm * x)
{
    vect3D temp;
    temp.absi = X.absi;
    temp.ordo = X.ordo;
    temp.cote = X.cote;
    x->u = rocv_scal_toqc(temp, U);
    x->v = rocv_scal_toqc(temp, V);
}


double wotf_stat_kefz(double lambda, point A1, point B1, point C1, point A2, point B2, point C2, int *idx)
{
    double res;
    parm a1, b1, c1, b2;
    vect3D N1, N2, N, U, V;
    qosp_unit_zamk(A1, B1, C1, &N1);
    qosp_unit_zamk(A2, B2, C2, &N2);
    res = rocv_scal_toqc(N1, N2);
    *idx = 2;
    if (res >= lambda) {
        N.absi = 0.5 * (N1.absi + N2.absi);
        N.ordo = 0.5 * (N1.ordo + N2.ordo);
        N.cote = 0.5 * (N1.cote + N2.cote);
        qubr_norm_foqk(&N);
        kowv_loca_colb(N, &U, &V);
        gufw_plan_mesv(U, V, A1, &a1);
        gufw_plan_mesv(U, V, B1, &b1);
        gufw_plan_mesv(U, V, C1, &c1);
        gufw_plan_mesv(U, V, B2, &b2);
        *idx = vefd_test_wacj(a1, b1, c1, b2);

    }
    return res;
}


void tukc_comm_nudj(manif_tl msh, int r, int s, int *a, int *b, int *forc_term)
{
    int i, j, k, n[3], m[3], *cmn;
    *forc_term = 0;
    n[0] = msh.entity[r].frvrt;
    n[1] = msh.entity[r].scvrt;
    n[2] = msh.entity[r].thvrt;
    m[0] = msh.entity[s].frvrt;
    m[1] = msh.entity[s].scvrt;
    m[2] = msh.entity[s].thvrt;
    cmn = (int *) malloc(9 * sizeof(int));
    k = 0;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            if (n[i] == m[j]) {
                cmn[k] = n[i];
                k++;
                break;
            }
        }
    if (k != 2) {
        fprintf(tmpout, "problem in adjacency\n");
        free(cmn);
        *forc_term = 1;
        return;
    }
    *a = cmn[0];
    *b = cmn[1];
    free(cmn);
}



double rocq_stat_todr(double lambda, manif_tl msh, int r, int s, int *idx, int *forc_term)
{
    int a1, b1, c1, a2, b2, c2, a, b, nx, f_trm;
    double res;
    *forc_term = 0;
    tukc_comm_nudj(msh, r, s, &a, &b, &f_trm);
    if (f_trm == 1) {
        *forc_term = 1;
        return 0.0;
    }

    nx = kuts_next_kolv(msh.entity[r], a);
    if (nx != b) {
        a1 = a;
        b1 = nx;
        c1 = b;
    } else {
        a1 = b;
        b1 = kuts_next_kolv(msh.entity[r], a1);
        c1 = a;
    }

    nx = kuts_next_kolv(msh.entity[s], a);
    if (nx != b) {
        c2 = b;
        a2 = a;
        b2 = kuts_next_kolv(msh.entity[s], a2);
    } else {
        c2 = a;
        a2 = b;
        b2 = kuts_next_kolv(msh.entity[s], a2);
    }

    res = wotf_stat_kefz(lambda, msh.knot[a1], msh.knot[b1], msh.knot[c1], msh.knot[a2], msh.knot[b2], msh.knot[c2], idx);
    return res;
}
