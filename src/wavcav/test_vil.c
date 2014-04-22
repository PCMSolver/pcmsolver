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
#include "geodesic.h"
#include "coarsequad.h"


void furk_prep_huqs(manif_tl * msh, int max_ned)
{
    fprintf(tmpout, "measure\n");
    hobf_fill_cewr(msh);
    fprintf(tmpout, "fill manifold\n");
    cogv_fill_zicd(msh, max_ned);
    mokq_chec_nukb(*msh);
    fprintf(tmpout, "manifold orientation\n");
    tujh_orie_novd(msh);
    fprintf(tmpout, "incidence\n");
    qosr_fill_fedt(msh);
}


void rudk_allo_tamq(int nnd, int nel, int max_ned, manif_tl * msh)
{
    msh->knot = (point *) malloc(nnd * sizeof(point));
    msh->entity = (telolf *) malloc(nel * sizeof(telolf));
    msh->increl = (teboka_topo *) malloc(nnd * sizeof(teboka_topo));
    msh->kt = (kt_t *) malloc(max_ned * sizeof(kt_t));
}


void lawn_dest_jukt(manif_tl * msh)
{
    free(msh->knot);
    free(msh->entity);
    free(msh->increl);
    free(msh->kt);
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

