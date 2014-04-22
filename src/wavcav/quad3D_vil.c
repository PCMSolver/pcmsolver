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
#include <stdio.h>
#include <stdlib.h>
#include "cavity.h"
#include "splinemol.h"
#include "pln_sph.h"
#include "triang.h"


int jipd_dete_getn(fajor_sion3D QUAD, kt_t * ed, int max_number)
{
    int i, nnd, nel, n1, n2, n3, n4, N, s, t, *nd_a, *nd_b, p;
    nd_a = (int *) malloc(4 * sizeof(int));
    nd_b = (int *) malloc(4 * sizeof(int));
    nnd = QUAD.n_grs;
    nel = QUAD.e_grs;
    for (i = 0; i < max_number; i++)
        ed[i].scent = -1;
    N = 0;
    for (i = 0; i < nel; i++) {
        n1 = QUAD.elem[i].frvrt;
        n2 = QUAD.elem[i].scvrt;
        n3 = QUAD.elem[i].thvrt;
        n4 = QUAD.elem[i].ftvrt;
        nd_a[0] = n1;
        nd_b[0] = n2;
        nd_a[1] = n2;
        nd_b[1] = n3;
        nd_a[2] = n3;
        nd_b[2] = n4;
        nd_a[3] = n4;
        nd_b[3] = n1;

        for (p = 0; p < 4; p++) {
            t = sujm_test_fujk(ed, nd_a[p], nd_b[p], N, &s);
            if (s == -1) {
                ed[N].frvrt = nd_a[p];
                ed[N].scvrt = nd_b[p];
                ed[N].frent = i;
                ed[N].scent = -1;
                N++;
            } else
                ed[t].scent = i;
        }
    }
    free(nd_a);
    free(nd_b);
    return N;
}


void teqr_fill_regm(fajor_sion3D * QUAD, int max_ned)
{
    int max_edge, nel, N, i, k;
    kt_t *ed;
    nel = QUAD->e_grs;
    max_edge = 4 * nel;
    ed = (kt_t *) malloc(max_edge * sizeof(kt_t));
    N = jipd_dete_getn(*QUAD, ed, max_edge);
    if (N >= max_ned) {
        fprintf(tmpout, "max_ned=%d is reached\n", max_ned);
        exit(0);
    }
    QUAD->k_grs = N;
    for (i = 0; i < N; i++) {
        QUAD->kt[i].frent = ed[i].frent;
        QUAD->kt[i].scent = ed[i].scent;
        QUAD->kt[i].frvrt = ed[i].frvrt;
        QUAD->kt[i].scvrt = ed[i].scvrt;
    }
    for (i = 0; i < QUAD->e_grs; i++) {
        QUAD->elem[i].frkt = -1;
        QUAD->elem[i].sckt = -1;
        QUAD->elem[i].trkt = -1;
        QUAD->elem[i].ftkt = -1;
    }
    for (i = 0; i < N; i++) {
        k = ed[i].frent;
        if (QUAD->elem[k].frkt == -1)
            QUAD->elem[k].frkt = i;
        else if (QUAD->elem[k].sckt == -1)
            QUAD->elem[k].sckt = i;
        else if (QUAD->elem[k].trkt == -1)
            QUAD->elem[k].trkt = i;
        else if (QUAD->elem[k].ftkt == -1)
            QUAD->elem[k].ftkt = i;
        k = ed[i].scent;
        if (k != -1) {
            if (QUAD->elem[k].frkt == -1)
                QUAD->elem[k].frkt = i;
            else if (QUAD->elem[k].sckt == -1)
                QUAD->elem[k].sckt = i;
            else if (QUAD->elem[k].trkt == -1)
                QUAD->elem[k].trkt = i;
            else if (QUAD->elem[k].ftkt == -1)
                QUAD->elem[k].ftkt = i;
        }
    }
    free(ed);
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

