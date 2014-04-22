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
#include <stdlib.h>
#include <stdlib.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "smooth.h"
#include "coarsequad.h"
#include "pln_sph.h"
#include "sas.h"
#include "meshsas.h"



void nirc_remo_tocp(manif_tl * msh, int *map_el, int nel_old)
{
    int nel, i, j, k, n1, n2, n3, ts1, ts2, ts3;
    int dummy, n_loc[3], ind, *map_loc, z;
    telolf *temp;
    nel = msh->e_grs;
    map_loc = (int *) malloc(nel * sizeof(int));
    temp = (telolf *) malloc(nel * sizeof(telolf));
    k = 0;
    for (i = 0; i < nel; i++) {
        n1 = msh->entity[i].frvrt;
        n2 = msh->entity[i].scvrt;
        n3 = msh->entity[i].thvrt;
        ind = 0;
        for (j = 0; j < k; j++) {
            n_loc[0] = temp[j].frvrt;
            n_loc[1] = temp[j].scvrt;
            n_loc[2] = temp[j].thvrt;
            ts1 = gonl_arra_govj(n_loc, 3, n1, &dummy);
            if (ts1 == 1) {
                ts2 = gonl_arra_govj(n_loc, 3, n2, &dummy);
                if (ts2 == 1) {
                    ts3 = gonl_arra_govj(n_loc, 3, n3, &dummy);
                    if (ts3 == 1) {
                        ind = 1;
                        break;
                    }
                }
            }
        }
        if (ind == 0) {
            temp[k].frvrt = n1;
            temp[k].scvrt = n2;
            temp[k].thvrt = n3;
            map_loc[i] = k;
            k++;
        }
    }

    for (i = 0; i < k; i++) {
        msh->entity[i].frvrt = temp[i].frvrt;
        msh->entity[i].scvrt = temp[i].scvrt;
        msh->entity[i].thvrt = temp[i].thvrt;
    }
    msh->e_grs = k;
    for (i = 0; i < nel_old; i++) {
        z = map_el[i];
        if (z != -1)
            map_el[i] = map_loc[z];
    }
    free(map_loc);
    free(temp);
}


void dopj_remo_qosd(manif_tl * msh)
{
    int nel, i, j, k, n1, n2, n3, ts1, ts2, ts3;
    int dummy, n_loc[3], ind;
    telolf *temp;
    nel = msh->e_grs;
    temp = (telolf *) malloc(nel * sizeof(telolf));
    k = 0;
    for (i = 0; i < nel; i++) {
        n1 = msh->entity[i].frvrt;
        n2 = msh->entity[i].scvrt;
        n3 = msh->entity[i].thvrt;
        ind = 0;
        for (j = 0; j < k; j++) {
            n_loc[0] = temp[j].frvrt;
            n_loc[1] = temp[j].scvrt;
            n_loc[2] = temp[j].thvrt;
            ts1 = gonl_arra_govj(n_loc, 3, n1, &dummy);
            if (ts1 == 1) {
                ts2 = gonl_arra_govj(n_loc, 3, n2, &dummy);
                if (ts2 == 1) {
                    ts3 = gonl_arra_govj(n_loc, 3, n3, &dummy);
                    if (ts3 == 1) {
                        ind = 1;
                        break;
                    }
                }
            }
        }
        if (ind == 0) {
            temp[k].frvrt = n1;
            temp[k].scvrt = n2;
            temp[k].thvrt = n3;
            k++;
        }
    }

    for (i = 0; i < k; i++) {
        msh->entity[i].frvrt = temp[i].frvrt;
        msh->entity[i].scvrt = temp[i].scvrt;
        msh->entity[i].thvrt = temp[i].thvrt;
    }
    msh->e_grs = k;
    free(temp);
    if (k != nel)
        fprintf(tmpout, "Nb dupl tri=%d\n", nel - k);
}


int sdlo = 0;



int fucm_disc_saqh(double scl, double scl2, manif_tl * msh, int max_ned, int *map_el)
{
    int *ed_col, *map_nd, nb_col, n_inv, *inv_node;
    int nnd_old, ned_old, max_inv, ned, suc = FAILURE;
    int max_col_aux = 100, max_col, nel_old, some[2] = { 802, 808 };
    nel_old = msh->e_grs;
    ned = msh->k_grs;
    if (ned < max_col_aux)
        max_col = ned;
    else
        max_col = max_col_aux;
    ed_col = (int *) malloc(max_col * sizeof(int));
    max_inv = 50 * max_col;
    if (max_inv > msh->n_grs)
        max_inv = msh->n_grs;
    ned = msh->k_grs;
    inv_node = (int *) malloc(max_inv * sizeof(int));
    nb_col = fahg_edge_lujd(scl, scl2, *msh, ed_col, max_col, inv_node, &n_inv, max_inv);
    fprintf(tmpout, "sdlo=%d    nb_col=%d\n", sdlo, nb_col);


    if (nb_col >= 1) {
        map_nd = (int *) malloc(msh->n_grs * sizeof(int));
        nnd_old = msh->n_grs;
        ned_old = msh->k_grs;

        gewn_chec_qamk(*msh);

        wacl_coll_hijt(msh, ed_col, nb_col, map_nd, map_el);


        cogv_fill_zicd(msh, max_ned);
        mokq_chec_nukb(*msh);
        qosr_fill_fedt(msh);

        gewn_chec_qamk(*msh);
        free(map_nd);
        suc = SUCCESS;
    }
    free(inv_node);
    free(ed_col);


    sdlo++;
    return suc;
}



void nusp_disc_jidk(double scl, double scl2, manif_tl * msh, rgb_lk * col, sphere * S, int max_ned)
{
    int nel, *map_el, suc, nb_trav = 100, i, j, nel_old, z;
    rgb_lk *temp;
    sphere *S_temp;
    nel = msh->e_grs;
    map_el = (int *) malloc((nel + 10) * sizeof(int));
    temp = (rgb_lk *) malloc((nel + 10) * sizeof(rgb_lk));
    S_temp = (sphere *) malloc((nel + 10) * sizeof(sphere));
    for (i = 0; i < nb_trav; i++) {
        fprintf(tmpout, "REMOVE BAD ELEMENTS(1)  TRAV=%d/%d\n", i, nb_trav - 1);
        nel_old = msh->e_grs;
        suc = fucm_disc_saqh(scl, scl2, msh, max_ned, map_el);
        if (suc == FAILURE) {
            fprintf(tmpout, "breaking\n");
            break;
        }
        if (suc == SUCCESS) {
            for (j = 0; j < nel_old; j++) {
                z = map_el[j];
                if (z != -1) {
                    temp[z].red = col[j].red;
                    temp[z].green = col[j].green;
                    temp[z].blue = col[j].blue;
                    neqg_find_lodr_bogm(S[j], &S_temp[z]);
                }
            }
            for (j = 0; j < msh->e_grs; j++) {
                col[j].red = temp[j].red;
                col[j].green = temp[j].green;
                col[j].blue = temp[j].blue;
                neqg_find_lodr_bogm(S_temp[j], &S[j]);
            }
        }
        if ((i == nb_trav - 1) && (suc == SUCCESS)) {
            fprintf(tmpout, "WARNING: nb_trav is reached but still some lavmah\n");

        }
        fprintf(tmpout, "-------------\n");
    }
    free(temp);
    free(map_el);
    free(S_temp);

}


void wihp_disc_sogj(double scl, double scl2, manif_tl * msh, int max_ned)
{
    int nel, *map_el, suc, nb_trav = 100, i, nel_old;
    nel = msh->e_grs;
    map_el = (int *) malloc((nel + 10) * sizeof(int));
    for (i = 0; i < nb_trav; i++) {
        fprintf(tmpout, "REMOVE BAD ELEMENTS(2)  TRAV=%d/%d\n", i, nb_trav - 1);
        nel_old = msh->e_grs;
        suc = fucm_disc_saqh(scl, scl2, msh, max_ned, map_el);
        if (suc == FAILURE) {
            fprintf(tmpout, "breaking\n");
            break;
        }
        if ((i == nb_trav - 1) && (suc == SUCCESS)) {
            fprintf(tmpout, "2-nb_trav is reached but still some lavmah\n");

        }
        mokq_chec_nukb(*msh);
        fprintf(tmpout, "-------------\n");
    }
    free(map_el);
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

