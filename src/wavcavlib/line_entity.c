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
#include <math.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"



void kivj_eval_gecz(line_entity le,double t,point *sol)
{sol->absi=le.p1.absi+t*(le.p2.absi-le.p1.absi);
sol->ordo=le.p1.ordo+t*(le.p2.ordo-le.p1.ordo);
sol->cote=le.p1.cote+t*(le.p2.cote-le.p1.cote);
}


void fagj_eval_sodh(line_entity le,double t,parm *sol)
{sol->u=le.p1.absi+t*(le.p2.absi-le.p1.absi);
sol->v=le.p1.ordo+t*(le.p2.ordo-le.p1.ordo);
}



void fenq_inve_dusj(line_entity C_in,line_entity *C_out)
{getf_find_rogc_todj(C_in.p1,&C_out->p2);
getf_find_rogc_todj(C_in.p2,&C_out->p1);
}


void cerv_disp_nods(line_entity le)
{fprintf(tmpout,"\tP1=(%f,%f,%f)\n",le.p1.absi,le.p1.ordo,le.p1.cote);
fprintf(tmpout,"\tP2=(%f,%f,%f)\n",le.p2.absi,le.p2.ordo,le.p2.cote);
}



 
