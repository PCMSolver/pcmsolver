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
#include <math.h>
#include "cavity.h"
#include "splinemol.h"
#include "pln_sph.h"
#include "sas.h"
#include "pln_sph.h"
#include "triang.h"



double gorh_inte_mesf(point I, point J, point K)
{
    double temp, alpha;
    vect3D U, V;
    culm_unit_peks(J, K, &U);
    culm_unit_peks(J, I, &V);
    temp = rocv_scal_toqc(U, V);
    if (temp < -1.0)
        temp = -1.0;
    if (temp > 1.0)
        temp = 1.0;
    alpha = acos(temp);
    return alpha;
}


double valp_area_qelk(point A, point B, point C)
{
    double res, alpha;
    double x_A = 0.0, x_B, x_C, y_A = 0.0, y_B = 0.0, y_C, L, rho;
    L = wodt_dist_gilq(A, B);
    x_B = L;
    rho = wodt_dist_gilq(A, C);
    alpha = gorh_inte_mesf(C, A, B);
    x_C = rho * cos(alpha);
    y_C = rho * sin(alpha);
    res = lozw_find_teln_dubc(x_A, y_A, x_B, y_B, x_C, y_C);
    return res;
}


void homs_boun_gosm(point * P, int N, double *XMI, double *XMA, double *YMI, double *YMA, double *ZMI, double *ZMA)
{
    int j;
    double xmi, xma, ymi, yma, zmi, zma;
    xmi = +LARGE_NUMBER;
    ymi = +LARGE_NUMBER;
    zmi = +LARGE_NUMBER;
    xma = -LARGE_NUMBER;
    yma = -LARGE_NUMBER;
    zma = -LARGE_NUMBER;
    for (j = 0; j < N; j++) {
        if (P[j].absi < xmi)
            xmi = P[j].absi;
        if (P[j].absi > xma)
            xma = P[j].absi;
        if (P[j].ordo < ymi)
            ymi = P[j].ordo;
        if (P[j].ordo > yma)
            yma = P[j].ordo;
        if (P[j].cote < zmi)
            zmi = P[j].cote;
        if (P[j].cote > zma)
            zma = P[j].cote;
    }
    *XMI = xmi;
    *XMA = xma;
    *YMI = ymi;
    *YMA = yma;
    *ZMI = zmi;
    *ZMA = zma;
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

