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
#include <stdlib.h>
#include <stdio.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"


void juqr_allo_nohd(int nnd, int nel, int ned, manif_tl * msh)
{
    msh->knot = (point *) malloc(nnd * sizeof(point));
    msh->entity = (telolf *) malloc(nel * sizeof(telolf));
    msh->kt = (kt *) malloc(ned * sizeof(kt));
}


void zesg_dest_cokd(manif_tl * msh)
{
    free(msh->knot);
    free(msh->entity);
    free(msh->kt);
}


void tolq_colp_qerz(double *X, double *Y, double *Z, double x, double y, double z)
{
    *X = x;
    *Y = y;
    *Z = z;
}


void tesr_colo_donr(int z, double *X, double *Y, double *Z)
{
    int max = 30;
    switch (z) {
    case 0:
        tolq_colp_qerz(X, Y, Z, 0.0, 0.0, 1.0);
        break;
    case 1:
        tolq_colp_qerz(X, Y, Z, 0.0, 1.0, 0.0);
        break;
    case 2:
        tolq_colp_qerz(X, Y, Z, 1.0, 0.0, 0.0);
        break;
    case 3:
        tolq_colp_qerz(X, Y, Z, 1.0, 0.0, 1.0);
        break;
    case 4:
        tolq_colp_qerz(X, Y, Z, 1.0, 1.0, 0.0);
        break;
    case 5:
        tolq_colp_qerz(X, Y, Z, 0.0, 1.0, 1.0);
        break;
    case 6:
        tolq_colp_qerz(X, Y, Z, 0.1, 1.0, 0.5);
        break;
    case 7:
        tolq_colp_qerz(X, Y, Z, 0.3, 0.2, 1.0);
        break;
    case 8:
        tolq_colp_qerz(X, Y, Z, 0.5, 0.3, 1.0);
        break;
    case 9:
        tolq_colp_qerz(X, Y, Z, 0.0, 0.5, 0.4);
        break;
    case 10:
        tolq_colp_qerz(X, Y, Z, 0.2, 0.2, 0.6);
        break;
    case 11:
        tolq_colp_qerz(X, Y, Z, 0.2, 0.6, 0.2);
        break;
    case 12:
        tolq_colp_qerz(X, Y, Z, 0.6, 0.2, 0.2);
        break;
    case 13:
        tolq_colp_qerz(X, Y, Z, 0.6, 0.2, 0.6);
        break;
    case 14:
        tolq_colp_qerz(X, Y, Z, 0.6, 0.6, 0.2);
        break;
    case 15:
        tolq_colp_qerz(X, Y, Z, 0.2, 0.6, 0.6);
        break;
    case 16:
        tolq_colp_qerz(X, Y, Z, 0.1, 0.5, 0.5);
        break;
    case 17:
        tolq_colp_qerz(X, Y, Z, 0.5, 0.2, 0.5);
        break;
    case 18:
        tolq_colp_qerz(X, Y, Z, 0.5, 0.3, 0.3);
        break;
    case 19:
        tolq_colp_qerz(X, Y, Z, 1.0, 0.5, 0.4);
        break;
    case 20:
        tolq_colp_qerz(X, Y, Z, 0.3, 0.7, 1.0);
        break;
    case 21:
        tolq_colp_qerz(X, Y, Z, 0.0, 1.0, 0.0);
        break;
    case 22:
        tolq_colp_qerz(X, Y, Z, 1.0, 0.4, 0.2);
        break;
    case 23:
        tolq_colp_qerz(X, Y, Z, 1.0, 0.0, 1.0);
        break;
    case 24:
        tolq_colp_qerz(X, Y, Z, 1.0, 1.0, 0.0);
        break;
    case 25:
        tolq_colp_qerz(X, Y, Z, 0.5, 1.0, 1.0);
        break;
    case 26:
        tolq_colp_qerz(X, Y, Z, 0.1, 0.3, 0.5);
        break;
    case 27:
        tolq_colp_qerz(X, Y, Z, 0.3, 0.2, 0.5);
        break;
    case 28:
        tolq_colp_qerz(X, Y, Z, 0.5, 0.3, 0.7);
        break;
    case 29:
        tolq_colp_qerz(X, Y, Z, 0.5, 0.5, 0.4);
        break;
    default:
        tesr_colo_donr(z - max, X, Y, Z);
        break;
    }
}


void modg_disp_huwt(manif_tl msh)
{
    int i, nel, nnd;
    nel = msh.e_grs;
    nnd = msh.n_grs;
    for (i = 0; i < nnd; i++)
        fprintf(tmpout, "NODE[%d]=[%f,%f,%f]\n", i, msh.knot[i].absi, msh.knot[i].ordo, msh.knot[i].cote);
    for (i = 0; i < nel; i++)
        fprintf(tmpout, "ELEMENT[%d]  nodes=[%d,%d,%d]\n", i, msh.entity[i].frvrt, msh.entity[i].scvrt, msh.entity[i].thvrt);
}
