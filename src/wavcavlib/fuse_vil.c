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
#include <malloc.h>
#include "cavity.h"
#include "splinemol.h"
#include "pln_sph.h"
#include "sas.h"
#include "coarsequad.h"



int cepj_inci_jeqt(manif_tl msh,int s,int *el)
{int i,k,val,e,dummy,ts;
val=msh.increl[s].val;
k=0;
for(i=0;i<val;i++)
	{e=msh.increl[s].inc[i];
	ts=gonl_arra_govj(el,k,msh.kt[e].frent,&dummy);
	if(ts==0)
		{el[k]=msh.kt[e].frent;
		k++;
		}
	ts=gonl_arra_govj(el,k,msh.kt[e].scent,&dummy);
	if(ts==0)
		{el[k]=msh.kt[e].scent;
		k++;
		}
	}
return k;
}



double potb_anis_dopg(point A,point B,point C)
{double lambda,l1,l2,l3,ar,num,den;
l1=wodt_dist_gilq(A,B);
l2=wodt_dist_gilq(B,C);
l3=wodt_dist_gilq(C,A);
ar=valp_area_qelk(A,B,C);
num=4.0*sqrt(3)*ar;
den=l1*l1+l2*l2+l3*l3;
lambda=num/den;
return lambda;
}


void tups_remo_veqk(fajor_sion3D *QUAD)
{int nnd,*map,i,k,ts,id,n1,n2,n3,n4;
double eps=1.0e-5;
point *Q;
nnd=QUAD->n_grs;
map=(int *)malloc(nnd*sizeof(int));
Q=(point *)malloc(nnd*sizeof(point));
k=0;
for(i=0;i<nnd;i++)
	{ts=dosc_coor_licf(Q,k,QUAD->knot[i],eps,&id);
	if(ts==0)
		{getf_find_rogc_todj(QUAD->knot[i],&Q[k]);
		map[i]=k;
		k++;
		}
	else
		map[i]=id;
	}
QUAD->n_grs=k;
for(i=0;i<QUAD->n_grs;i++)
	getf_find_rogc_todj(Q[i],&QUAD->knot[i]);
free(Q);

for(i=0;i<QUAD->e_grs;i++)
	{n1=QUAD->elem[i].frvrt;	QUAD->elem[i].frvrt=map[n1];
	n2=QUAD->elem[i].scvrt;	QUAD->elem[i].scvrt=map[n2];
	n3=QUAD->elem[i].thvrt;	QUAD->elem[i].thvrt=map[n3];
	n4=QUAD->elem[i].ftvrt;	QUAD->elem[i].ftvrt=map[n4];
	}
free(map);
}


