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


void furk_prep_huqs(manif_tl *msh,int max_ned)
{fprintf(tmpout,"measure\n");
hobf_fill_cewr(msh);
fprintf(tmpout,"fill manifold\n");
cogv_fill_zicd(msh,max_ned);
mokq_chec_nukb(*msh);
fprintf(tmpout,"manifold orientation\n");
tujh_orie_novd(msh);
fprintf(tmpout,"incidence\n");
qosr_fill_fedt(msh);
}


void rudk_allo_tamq(int nnd,int nel,int max_ned,manif_tl *msh)
{msh->knot=(point *)malloc(nnd*sizeof(point));
msh->entity=(telolf *)malloc(nel*sizeof(telolf));
msh->increl=(teboka_topo *)malloc(nnd*sizeof(teboka_topo));
msh->kt=(kt *)malloc(max_ned*sizeof(kt));
}

 
void lawn_dest_jukt(manif_tl *msh)
{free(msh->knot);
free(msh->entity);
free(msh->increl);
free(msh->kt);
}

