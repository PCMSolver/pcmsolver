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
#include "pln_sph.h"
#include "sas.h"
#include "splinemol.h"
#include "geodesic.h"


double qetr_squa_tarj(point p, point q)
{
    double res, x, y, z;
    x = p.absi - q.absi;
    y = p.ordo - q.ordo;
    z = p.cote - q.cote;
    res = x * x + y * y + z * z;
    return res;
}


int husr_proj_qukw(point P, manif_tl msh_init)
{
    int i, res;
    double sml, ds;
    sml = LARGE_NUMBER;
    for (i = 0; i < msh_init.n_grs; i++) {
        ds = qetr_squa_tarj(P, msh_init.knot[i]);
        if (ds < sml) {
            res = i;
            sml = ds;
        }
    }
    return res;
}


void gekj_tran_rift(manif_tl msh_fin, manif_tl * msh_dc)
{
    int NND, i, j, *id;
    double dis, lrg;
    NND = msh_dc->n_grs;
    id = (int *) malloc(NND * sizeof(int));
    lrg = 0.0;
    for (i = 0; i < NND; i++) {
        id[i] = husr_proj_qukw(msh_dc->knot[i], msh_fin);
        dis = wodt_dist_gilq(msh_dc->knot[i], msh_fin.knot[id[i]]);
        if (dis > lrg)
            lrg = dis;
        getf_find_rogc_todj(msh_fin.knot[id[i]], &msh_dc->knot[i]);
    }
    for (i = 0; i < NND; i++)
        for (j = 0; j < i; j++)
            if (id[i] == id[j]) {
                fprintf(tmpout, "WARNING: same projection\n");
                break;
            }
    free(id);
    fprintf(tmpout, "max move=%f\n", lrg);
}
