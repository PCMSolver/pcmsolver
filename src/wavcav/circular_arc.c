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

#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
//#include  <malloc.h>
#include "cavity.h"
#include  "pln_sph.h"
#include  "eval.h"



void vewd_comp_tilg(c_arc ca, double *t2, double *t3)
{
    double R, temp1, temp2, cs, sn, pi, T2, T3;
    temp1 = ca.p2.absi - ca.p1.absi;
    temp2 = ca.p2.ordo - ca.p1.ordo;
    R = sqrt(temp1 * temp1 + temp2 * temp2);
    cs = (ca.p2.absi - ca.p1.absi) / R;
    sn = (ca.p2.ordo - ca.p1.ordo) / R;
    pi = 3.14159265358979323846264338328;
    if (sn >= 0.0)
        T2 = acos(cs);
    else
        T2 = 2.0 * pi - acos(cs);

    cs = (ca.p3.absi - ca.p1.absi) / R;
    sn = (ca.p3.ordo - ca.p1.ordo) / R;
    if (sn >= 0.0)
        T3 = acos(cs);
    else
        T3 = 2.0 * pi - acos(cs);
    *t2 = T2;
    if (T3 >= T2)
        *t3 = T3;
    else
        *t3 = T3 + 2 * pi;
}



void qupl_find_bofh_honw(double **A, double *b, double *c, int n, int m)
{
    int i, j;
    double temp;
    for (i = 0; i < n; i++) {
        temp = 0.0;
        for (j = 0; j < m; j++)
            temp = temp + A[i][j] * b[j];
        c[i] = temp;
    }
}


void polt_appl_remc(xform OP, point x, point * y)
{
    int i, j;
    double *temp, *out, **MAT;
    temp = (double *) malloc(3 * sizeof(double));
    out = (double *) malloc(3 * sizeof(double));
    temp[0] = x.absi;
    temp[1] = x.ordo;
    temp[2] = x.cote;
    MAT = allocate_mat(3, 3);
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            MAT[i][j] = OP.R[i][j];
    qupl_find_bofh_honw(MAT, temp, out, 3, 3);
    y->absi = out[0] + OP.T[0];
    y->ordo = out[1] + OP.T[1];
    y->cote = out[2] + OP.T[2];
    free(temp);
    free(out);
    tehg_free_dacp(MAT, 3, 3);
}


void bojv_appl_locq(xform OP, point x, parm * y)
{
    int i, j;
    double *temp, *out, **MAT;
    temp = (double *) malloc(3 * sizeof(double));
    out = (double *) malloc(2 * sizeof(double));
    temp[0] = x.absi;
    temp[1] = x.ordo;
    temp[2] = x.cote;
    MAT = allocate_mat(2, 3);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 3; j++)
            MAT[i][j] = OP.R[i][j];
    qupl_find_bofh_honw(MAT, temp, out, 2, 3);
    y->u = out[0] + OP.T[0];
    y->v = out[1] + OP.T[1];
    free(temp);
    free(out);
    tehg_free_dacp(MAT, 2, 3);
}


void cowm_oppo_retq(xform OP_in, xform * OP_out)
{
    int i, j, k;
    double SYM[3][3];
    for (i = 0; i < 3; i++)
        OP_out->T[i] = OP_in.T[i];
    SYM[0][0] = 1.0;
    SYM[0][1] = 0.0;
    SYM[0][2] = 0.0;
    SYM[1][0] = 0.0;
    SYM[1][1] = -1.0;
    SYM[1][2] = 0.0;
    SYM[2][0] = 0.0;
    SYM[2][1] = 0.0;
    SYM[2][2] = 1.0;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            OP_out->R[i][j] = 0.0;
            for (k = 0; k < 3; k++)
                OP_out->R[i][j] = OP_out->R[i][j] + OP_in.R[i][k] + SYM[k][j];
        }
}


void petl_eval_rebm(c_arc ca, double t, point * sol)
{
    double R, temp1, temp2;
    point temp;
    temp1 = ca.p2.absi - ca.p1.absi;
    temp2 = ca.p2.ordo - ca.p1.ordo;
    R = sqrt(temp1 * temp1 + temp2 * temp2);
    temp.absi = ca.p1.absi + R * cos(t);
    temp.ordo = ca.p1.ordo + R * sin(t);
    temp.cote = ca.ZT;
    polt_appl_remc(ca.TR, temp, sol);
}



void pukc_eval_huvl(c_arc ca, double t, parm * sol)
{
    double R, temp1, temp2;
    point temp;
    temp1 = ca.p2.absi - ca.p1.absi;
    temp2 = ca.p2.ordo - ca.p1.ordo;
    R = sqrt(temp1 * temp1 + temp2 * temp2);
    temp.absi = ca.p1.absi + R * cos(t);
    temp.ordo = ca.p1.ordo + R * sin(t);
    temp.cote = ca.ZT;
    bojv_appl_locq(ca.TR, temp, sol);
}



void wunb_inve_zolr(c_arc C_in, c_arc * C_out)
{
    point start, term;
    getf_find_rogc_todj(C_in.p1, &C_out->p1);
    C_out->A = C_in.A;
    C_out->B = C_in.B;

    start.absi = C_in.p2.absi - C_in.p1.absi;
    start.ordo = C_in.p2.ordo - C_in.p1.ordo;
    start.cote = C_in.p2.cote - C_in.p1.cote;

    term.absi = C_in.p3.absi - C_in.p1.absi;
    term.ordo = C_in.p3.ordo - C_in.p1.ordo;
    term.cote = C_in.p3.cote - C_in.p1.cote;

    C_out->p2.absi = +term.absi + C_in.p1.absi;
    C_out->p2.ordo = -term.ordo + C_in.p1.ordo;
    C_out->p2.cote = +term.cote + C_in.p1.cote;

    C_out->p3.absi = +start.absi + C_in.p1.absi;
    C_out->p3.ordo = -start.ordo + C_in.p1.ordo;
    C_out->p3.cote = +start.cote + C_in.p1.cote;

    cowm_oppo_retq(C_in.TR, &C_out->TR);
}


void jumn_disp_cugw(c_arc ca)
{
    int i, j;
    fprintf(tmpout, "\tCenter point=(%f,%f,%f)\n", ca.p1.absi, ca.p1.ordo, ca.p1.cote);
    fprintf(tmpout, "\tStart point=(%f,%f,%f)\n", ca.p2.absi, ca.p2.ordo, ca.p2.cote);
    fprintf(tmpout, "\tEnd point=(%f,%f,%f)\n", ca.p3.absi, ca.p3.ordo, ca.p3.cote);
    fprintf(tmpout, "\tz-value=%f\n", ca.ZT);
    fprintf(tmpout, "\toperator:\n");
    for (i = 0; i < 3; i++) {
        fprintf(tmpout, "\t\t");
        for (j = 0; j < 3; j++)
            fprintf(tmpout, "%f  ", ca.TR.R[i][j]);
        fprintf(tmpout, "%f\n", ca.TR.T[i]);
    }
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

