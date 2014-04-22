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
#include <math.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"
#include "splinemol.h"
#include "meshsas.h"
#include "splinemol.h"


void hobf_fill_cewr(manif_tl * msh)
{
    int nel, i, n1, n2, n3;
    nel = msh->e_grs;
    for (i = 0; i < nel; i++) {
        n1 = msh->entity[i].frvrt;
        n2 = msh->entity[i].scvrt;
        n3 = msh->entity[i].thvrt;
        msh->entity[i].ms_wert = valp_area_qelk(msh->knot[n1], msh->knot[n2], msh->knot[n3]);
    }
}



void qosr_fill_fedt(manif_tl * msh)
{
    int i, j, k, ts, nnd, ned, n[2], vl, m, mxvl;
    nnd = msh->n_grs;
    ned = msh->k_grs;
    for (i = 0; i < nnd; i++)
        msh->increl[i].val = 0;
    for (i = 0; i < ned; i++) {
        n[0] = msh->kt[i].frvrt;
        n[1] = msh->kt[i].scvrt;
        for (j = 0; j < 2; j++) {
            m = n[j];
            vl = msh->increl[m].val;
            ts = -1;
            for (k = 0; k < vl; k++)
                if (msh->increl[m].inc[k] == i) {
                    ts = 1;
                    break;
                }
            if (ts == -1) {
                if (vl >= MAXVALENCE) {
                    fprintf(tmpout, "MAXVALENCE=%d is reached\n", MAXVALENCE);
                    exit(0);
                }
                msh->increl[m].inc[vl] = i;
                msh->increl[m].val = vl + 1;
            }
        }
    }
    mxvl = 0;
    for (i = 0; i < msh->n_grs; i++)
        if (msh->increl[i].val > mxvl)
            mxvl = msh->increl[i].val;
}



void degr_find_fuqd(manif_tl msh, double *XMIN, double *XMAX, double *YMIN, double *YMAX, double *ZMIN, double *ZMAX)
{
    int i, nel, nnd;
    double x, y, z, xmi, xma, ymi, yma, zmi, zma, lrg = LARGE_NUMBER;
    nel = msh.e_grs;
    nnd = msh.n_grs;
    xmi = lrg;
    xma = -lrg;
    ymi = lrg;
    yma = -lrg;
    zmi = lrg;
    zma = -lrg;
    if ((nel != 0) && (nnd != 0)) {
        for (i = 0; i < nnd; i++) {
            x = msh.knot[i].absi;
            y = msh.knot[i].ordo;
            z = msh.knot[i].cote;
            if (x < xmi)
                xmi = x;
            if (y < ymi)
                ymi = y;
            if (z < zmi)
                zmi = z;
            if (x > xma)
                xma = x;
            if (y > yma)
                yma = y;
            if (z > zma)
                zma = z;
        }
    } else {
        xmi = 0.0;
        xma = 0.0;
        ymi = 0.0;
        yma = 0.0;
        zmi = 0.0;
        zma = 0.0;
    }
    *XMIN = xmi;
    *XMAX = xma;
    *YMIN = ymi;
    *YMAX = yma;
    *ZMIN = zmi;
    *ZMAX = zma;
}


void goft_find_doqk(manif_tl msh, double *XMIN, double *XMAX, double *YMIN, double *YMAX, double *ZMIN, double *ZMAX)
{
    double xmin, xmax, ymin, ymax, zmin, zmax, eps = 0.01, diff;
    degr_find_fuqd(msh, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);
    diff = xmax - xmin;
    if (diff < eps) {
        *XMAX = xmax + 0.5;
        *XMIN = xmin - 0.5;
    } else {
        *XMAX = xmax;
        *XMIN = xmin;
    }

    diff = ymax - ymin;
    if (diff < eps) {
        *YMAX = ymax + 0.5;
        *YMIN = ymin - 0.5;
    } else {
        *YMAX = ymax;
        *YMIN = ymin;
    }

    diff = zmax - zmin;
    if (diff < eps) {
        *ZMAX = zmax + 0.5;
        *ZMIN = zmin - 0.5;
    } else {
        *ZMAX = zmax;
        *ZMIN = zmin;
    }
}



void nahd_find_vukd_gaph(manif_tl * msh, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
    int i, nnd, nel, n1, n2, n3;
    double len[3], lrg, den;
    point G;
    nnd = msh->n_grs;
    G.absi = 0.5 * (xmin + xmax);
    G.ordo = 0.5 * (ymin + ymax);
    G.cote = 0.5 * (zmin + zmax);
    for (i = 0; i < nnd; i++) {
        msh->knot[i].absi = msh->knot[i].absi - G.absi;
        msh->knot[i].ordo = msh->knot[i].ordo - G.ordo;
        msh->knot[i].cote = msh->knot[i].cote - G.cote;
    }

    len[0] = xmax - xmin;
    len[1] = ymax - ymin;
    len[2] = zmax - zmin;
    lrg = 0.0;
    for (i = 0; i < 3; i++)
        if (len[i] >= lrg)
            lrg = len[i];
    den = 0.5 * lrg;
    for (i = 0; i < nnd; i++) {
        msh->knot[i].absi = msh->knot[i].absi / den;
        msh->knot[i].ordo = msh->knot[i].ordo / den;
        msh->knot[i].cote = msh->knot[i].cote / den;
    }

    nel = msh->e_grs;
    for (i = 0; i < nel; i++) {
        n1 = msh->entity[i].frvrt;
        n2 = msh->entity[i].scvrt;
        n3 = msh->entity[i].thvrt;
        msh->entity[i].ms_wert = valp_area_qelk(msh->knot[n1], msh->knot[n2], msh->knot[n3]);
    }
}



int vist_find_jetl(manif_tl msh, int el, int *bulk, hash_entry * H)
{
    int id = -1, nel, i, p;
    nel = msh.e_grs;
    for (p = 0; p < H[el].nb; p++) {
        i = H[el].list[p];
        if ((i < 0) || (i >= nel)) {
            fprintf(tmpout, "Access beyond limit: i=%d   el=%d\n", i, el);
            exit(0);
        }
        if (bulk[i] == +1) {
            id = i;
            break;
        }
    }
    return id;
}


void wipz_inci_gemk(manif_tl msh, hash_entry * H)
{
    int nel, i, j, e[3], E1, E2;
    nel = msh.e_grs;
    for (i = 0; i < nel; i++)
        H[i].nb = 3;
    for (i = 0; i < nel; i++) {
        e[0] = msh.entity[i].frkt;
        e[1] = msh.entity[i].sckt;
        e[2] = msh.entity[i].trkt;
        for (j = 0; j < 3; j++) {
            E1 = msh.kt[e[j]].frent;
            E2 = msh.kt[e[j]].scent;
            if ((E1 >= nel) || (E2 >= nel)) {
                fprintf(tmpout, "incidence too large\n");
                exit(0);
            }
            if (E1 == i)
                H[i].list[j] = E2;
            else if (E2 == i)
                H[i].list[j] = E1;
            else {
                fprintf(tmpout, "Neither cases\n");
                exit(0);
            }
        }
    }
}


void jacm_chec_lejt(manif_tl msh)
{
    int i, ned, el1, el2;
    ned = msh.k_grs;
    for (i = 0; i < ned; i++) {
        el1 = msh.kt[i].frent;
        el2 = msh.kt[i].scent;
        if ((el1 == -1) || (el2 == -1)) {
            fprintf(tmpout, "minus sign boundary in manifold\n");
            exit(0);
        }
    }
    fprintf(tmpout, "MANIFOLD IS WELL CLOSED\n");
}



void tujh_orie_novd(manif_tl * msh)
{
    int nel, ned, *bulk, i, p, suc, ort, e1, e2, seed = 0, n_bk;
    hash_entry *H;
    jacm_chec_lejt(*msh);
    nel = msh->e_grs;
    H = (hash_entry *) malloc(nel * sizeof(hash_entry));
    for (i = 0; i < nel; i++)
        H[i].list = (int *) malloc(3 * sizeof(int));
    fprintf(tmpout, "incidence look:  nel=%d\n", nel);
    wipz_inci_gemk(*msh, H);
    bulk = (int *) malloc(nel * sizeof(int));
    bulk[seed] = +1;
    for (i = 0; i < nel; i++)
        if (i != seed)
            bulk[i] = -1;
    n_bk = 1;
    ned = msh->k_grs;
    fprintf(tmpout, "expansion\n");
    for (p = 0; p < ned; p++) {
        suc = qunf_expa_foqg(msh, bulk, H);
        if (suc == FAILURE)
            break;
    }
    for (i = 0; i < nel; i++)
        if (bulk[i] == -1) {
            fprintf(tmpout, "Expansion is not complete\n");
            exit(0);
        }
    fprintf(tmpout, "Complete expansion\n");
    free(bulk);
    cogv_fill_zicd(msh, ned);
    for (i = 0; i < msh->k_grs; i++) {
        e1 = msh->kt[i].frent;
        e2 = msh->kt[i].scent;
        ort = kotc_test_kolj(*msh, e1, e2);
        if (ort == 0) {
            fprintf(tmpout, "kt=%d\n", i);
            fprintf(tmpout, "WARNING: Inconsistent orientation\n");
            fprintf(tmpout, "element[%d]=[%d,%d,%d]   ", e1, msh->entity[e1].frvrt, msh->entity[e1].scvrt, msh->entity[e1].thvrt);
            fprintf(tmpout, "element[%d]=[%d,%d,%d]\n", e2, msh->entity[e2].frvrt, msh->entity[e2].scvrt, msh->entity[e2].thvrt);
            exit(0);
        }
    }
    for (i = 0; i < nel; i++)
        free(H[i].list);
    free(H);
    fprintf(tmpout, "GOOD LOCAL ORIENTATIONS\n");
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

