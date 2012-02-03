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
#include <math.h>
#include <stdlib.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"



void wolf_eval_murg(trmsrf ts,
double u,double v,point *sol)
{int TP,bdr;
parm p;
TP=ts.type;
bdr=ts.boundary;
switch(TP)
	{case 3:
		dufj_eval_wejf(ts.ts,u,v,sol);
		break;
	case 4:
		p.u=u;
		p.v=v;
		rogv_eval_dukw(ts.pt,p,sol);
		break;
	case 5:
		p.u=u;
		p.v=v;
		rows_eval_qusg(ts.pc,p,sol);
		break;
	case 6:
		p.u=u;
		p.v=v;
		cutn_eval_mecn(ts.bn,p,sol);
		break;
	}
}
 

void weht_disp_tosb(pt_tor PT)
{fprintf(tmpout,"ALPHA:\n");
nepf_disp_bulp(PT.alpha);
fprintf(tmpout,"BETA:\n");
nepf_disp_bulp(PT.beta);
fprintf(tmpout,"GAMMA:\n");
nepf_disp_bulp(PT.gamma);
fprintf(tmpout,"DELTA:\n");
nepf_disp_bulp(PT.delta);
}

 
void hewt_disp_mohw(pt_cnv PC)
{fprintf(tmpout,"sitn=%d\n",PC.sitn);
fprintf(tmpout,"ALPHA:\n");
nepf_disp_bulp(PC.alpha);
fprintf(tmpout,"BETA:\n");
nepf_disp_bulp(PC.beta);
fprintf(tmpout,"GAMMA:\n");
nepf_disp_bulp(PC.gamma);
}


void matc_disp_huqm(trmsrf ts)
{int tp,bdr,nb,i;
tp=ts.type;
bdr=ts.boundary;
nb=ts.nb_inner;
fprintf(tmpout,"BASE SURFACE INFORMATION:\n");
switch(tp)
	{case 3:
		fprintf(tmpout,"\tBase surface=trimmed sphere\n\n");
		jehg_disp_fecm(ts.ts);
		break;
	case 4:
		fprintf(tmpout,"\tBase surface=toroidal patch\n\n");
		weht_disp_tosb(ts.pt);
		break;
	case 5:
		fprintf(tmpout,"\tBase surface=concave patch\n\n");
		hewt_disp_mohw(ts.pc);
		break;
	}
fprintf(tmpout,"OUTER BOUNDARY INFORMATION:\n");
fprintf(tmpout,"\tType of outer boundary=%d\n",bdr);
if(bdr==1)
  {hewr_disp_toqh(ts.cc);
  }
fprintf(tmpout,"\tNumber of inner boundaries=%d\n",nb);
for(i=0;i<nb;i++)
	{fprintf(tmpout,"%d-th INNER BOUNDARY INFORMATION:\n",i);
	hewr_disp_toqh(ts.inner[i]);
	}
}


