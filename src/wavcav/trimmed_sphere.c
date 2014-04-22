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
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"
#include "eval.h"



void degj_para_mojc(trm_sph S, double u, double v, point * P)
{
    double den;
    den = 1.0 + u * u + v * v;
    P->absi = S.rad * 2.0 * u / den;
    P->ordo = S.rad * 2.0 * v / den;
    P->cote = S.rad * (1.0 - u * u - v * v) / den;
}



void peht_find_ligd_mebn(vect3D W_in, double phi, double theta, vect3D * W_out)
{
    int i;
    double **A, *V_in, *V_out;
    vect3D temp;

    A = (double **) malloc(3 * sizeof(double *));
    for (i = 0; i < 3; i++)
        A[i] = (double *) malloc(3 * sizeof(double));
    A[0][0] = 0.0;
    A[0][1] = 0.0;
    A[0][2] = 1.0;
    A[1][0] = 0.0;
    A[1][1] = 1.0;
    A[1][2] = 0.0;
    A[2][0] = -1.0;
    A[2][1] = 0.0;
    A[2][2] = 0.0;
    V_in = (double *) malloc(3 * sizeof(double));
    V_out = (double *) malloc(3 * sizeof(double));
    V_in[0] = W_in.absi;
    V_in[1] = W_in.ordo;
    V_in[2] = W_in.cote;
    fisd_find_rumj_veft(A, V_in, V_out);
    temp.absi = V_out[0];
    temp.ordo = V_out[1];
    temp.cote = V_out[2];

    mopb_dire_woqp(temp, phi, theta, W_out);
    for (i = 0; i < 3; i++)
        free(A[i]);
    free(A);
    free(V_in);
    free(V_out);
}


void sawk_find_vugj_mels(double phi, double theta, vect3D A, vect3D B)
{
    double d;
    vect3D temp;
    peht_find_ligd_mebn(B, phi, theta, &temp);
    fprintf(tmpout, "z_A=[%f,%f,%f]\n", A.absi, A.ordo, A.cote);
    fprintf(tmpout, "z_temp=[%f,%f,%f]\n", temp.absi, temp.ordo, temp.cote);
    d = wodt_dist_gilq(A, temp);
    if (d > 0.01) {
        fprintf(tmpout, "z_d=%f\n", d);
        exit(0);
    }
}



void tulr_find_rads_tojm(vect3D S_in, double phi, double theta, vect3D * S_out)
{
    int i;
    double **A, *V_in, *V_out;
    vect3D tp;
    V_in = (double *) malloc(3 * sizeof(double));
    V_out = (double *) malloc(3 * sizeof(double));
    fesg_inve_pahj(S_in, phi, theta, &tp);
    V_out[0] = tp.absi;
    V_out[1] = tp.ordo;
    V_out[2] = tp.cote;

    A = (double **) malloc(3 * sizeof(double *));
    for (i = 0; i < 3; i++)
        A[i] = (double *) malloc(3 * sizeof(double));
    A[0][0] = 0.0;
    A[0][1] = 0.0;
    A[0][2] = -1.0;
    A[1][0] = 0.0;
    A[1][1] = 1.0;
    A[1][2] = 0.0;
    A[2][0] = 1.0;
    A[2][1] = 0.0;
    A[2][2] = 0.0;
    fisd_find_rumj_veft(A, V_out, V_in);
    S_out->absi = V_in[0];
    S_out->ordo = V_in[1];
    S_out->cote = V_in[2];
    for (i = 0; i < 3; i++)
        free(A[i]);
    free(A);
    free(V_in);
    free(V_out);

}



void dufj_eval_wejf(trm_sph S, double u, double v, point * P)
{
    double phi, theta;
    point temp;
    vect3D h;

    degj_para_mojc(S, u, v, &temp);
    vewr_sphe_ruhd(S.nrml.absi, S.nrml.ordo, S.nrml.cote, &phi, &theta);
    peht_find_ligd_mebn(temp, phi, theta, &h);

    P->absi = h.absi + S.zent.absi;
    P->ordo = h.ordo + S.zent.ordo;
    P->cote = h.cote + S.zent.cote;
}


void jehg_disp_fecm(trm_sph T)
{
    fprintf(tmpout, "Trimmed sphere:\n");
    fprintf(tmpout, "\tcenter=[%f,%f,%f]\n", T.zent.absi, T.zent.ordo, T.zent.cote);
    fprintf(tmpout, "\tradius=%f\n", T.rad);
    fprintf(tmpout, "\tbeta=%f\n", T.beta);
    fprintf(tmpout, "\tnormal=[%f,%f,%f]\n", T.nrml.absi, T.nrml.ordo, T.nrml.cote);
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

