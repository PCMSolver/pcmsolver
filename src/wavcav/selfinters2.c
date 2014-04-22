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
#include "pln_sph.h"
#include "eval.h"
#include "sas.h"


void gapl_allo_necf(int nnd, int nel, conc_model * CM)
{
    CM->msh.knot = (point *) malloc(nnd * sizeof(point));
    CM->msh.entity = (telolf *) malloc(nel * sizeof(telolf));
    CM->G = (point *) malloc(nel * sizeof(point));
    CM->nrm = (vect3D *) malloc(nel * sizeof(vect3D));
}


void wabf_dest_purl(conc_model * CM)
{
    free(CM->msh.knot);
    free(CM->msh.entity);
    free(CM->G);
    free(CM->nrm);
}


void volq_find_luvc_jurq(trmsrf * surf3, int nb_surf3)
{
    int nnd, nel, i, j, k, *exc, p, q, N = 7;
    int ts_interf, nb_pairs, *some;
    int id_0, id_1, id_2;
    kt_t *E;
    manif_ro msh_tri;
    conc_model *CM;

    nnd = 3 * N * N;
    nel = 6 * (N - 1) * (N - 1);
    CM = (conc_model *) malloc(nb_surf3 * sizeof(conc_model));
    for (i = 0; i < nb_surf3; i++)
        gapl_allo_necf(nnd, nel, &CM[i]);
    msh_tri.knot = (parm *) malloc(nnd * sizeof(parm));
    msh_tri.entity = (telolf *) malloc(nel * sizeof(telolf));
    fiwh_simm_tucs(N, &msh_tri, &id_0, &id_1, &id_2);
    for (i = 0; i < nb_surf3; i++)
        kotr_find_tujg(msh_tri, id_0, id_1, id_2, surf3[i].pc, &CM[i]);

    free(msh_tri.knot);
    free(msh_tri.entity);
    fprintf(tmpout, "All triangulations are found\n");

    exc = (int *) malloc(nb_surf3 * sizeof(int));
    E = (kt_t *) malloc(nb_surf3 * sizeof(kt_t));
    for (i = 0; i < nb_surf3; i++)
        exc[i] = 0;
    k = 0;
    for (i = 0; i < nb_surf3; i++)
        if (exc[i] == 0) {
            for (j = 0; j < i; j++)
                if (exc[j] == 0) {
                    ts_interf = ticj_inte_ludj(CM[i], CM[j]);
                    if (ts_interf == 1) {
                        fprintf(tmpout, "Detected self-intersection: patches[%d,%d]\n", i, j);

                        E[k].frvrt = i;
                        E[k].scvrt = j;
                        k++;
                        exc[i] = 1;
                        exc[j] = 1;
                    }
                }
        }
    nb_pairs = k;
    free(exc);
    fprintf(tmpout, "Number of self-interfersections=%d\n", nb_pairs);

    some = (int *) malloc(2 * sizeof(int));
    for (i = 0; i < nb_pairs; i++) {
        p = E[i].frvrt;
        q = E[i].scvrt;
        surf3[p].pc.sitn = 1;
        surf3[q].pc.sitn = 1;
        surf3[p].pc.scl = 10.5;
        surf3[q].pc.scl = 10.5;

    }
    free(some);

    free(E);
    for (i = 0; i < nb_surf3; i++)
        wabf_dest_purl(&CM[i]);
    free(CM);
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

