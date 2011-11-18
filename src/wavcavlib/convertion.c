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


int lezc_node_bils(manif_tl msh,int s,int nd)
{int res;
res=0;
if(msh.entity[s].frvrt==nd)		res=1;
if(msh.entity[s].scvrt==nd)		res=1;
if(msh.entity[s].thvrt==nd)		res=1;
return res;
}



int falq_trav_noth(manif_tl msh,int e,int nd)
{int el[2],res=-1,i,ts;
el[0]=msh.kt[e].frent;
el[1]=msh.kt[e].scent;
for(i=0;i<2;i++)
	{ts=lezc_node_bils(msh,el[i],nd);
	if(ts==1)
		{res=el[i];
		break;
		}
	}
return res;
}



int mifz_comm_necg(manif_tl msh,int e,int f)
{int res=-1,E[2],F[2],i,j;
E[0]=msh.kt[e].frent;
E[1]=msh.kt[e].scent;
F[0]=msh.kt[f].frent;
F[1]=msh.kt[f].scent;
for(i=0;i<2;i++)
for(j=0;j<2;j++)
if(E[i]==F[j])
	{res=E[i];
	break;
	}
return res;
}



int gewp_supp_reqs(prat_main GR,int na,int nb)
{int ed=-1,val,i,e,n1,n2,nx;
val=GR.dgr[na];
for(i=0;i<val;i++)
	{e=GR.incd[na][i];
	n1=GR.kt[e].str;
	n2=GR.kt[e].ter;
	if(n1==na)	nx=n2;
	if(n2==na)	nx=n1;
	if(nx==nb)
		{ed=e;
		break;
		}
	}
return ed;
}



void with_conv_davf(manif_tl msh,prat_main GR,point *omega,
int *type,int *underl,int nb,float_curve *F)
{int i,tp_cr,tp_nx,e,tau,c_inc;

for(i=0;i<nb;i++)
	{F->stn[i].absi=omega[i].absi;
	F->stn[i].ordo=omega[i].ordo;
	F->stn[i].cote=omega[i].cote;
	if(type[i]==1)		F->kt_idx[i]=underl[i];
	if(type[i]==2)		F->nd_idx[i]=underl[i];
	}

for(i=0;i<nb-1;i++)
	{tp_cr=type[i];
	tp_nx=type[i+1];
	if((tp_cr==2)&&(tp_nx==2))		F->cs[i]=1;
	if((tp_cr==2)&&(tp_nx==1))		F->cs[i]=2;
	if((tp_cr==1)&&(tp_nx==1))		F->cs[i]=3;
	if((tp_cr==1)&&(tp_nx==2))		F->cs[i]=4;
	}

for(i=0;i<nb-1;i++)
	{if(F->cs[i]==1)
		{e=gewp_supp_reqs(GR,underl[i],underl[i+1]);
		if(e==-1)
			{fprintf(tmpout,"Unable to find supporing edge\n");
			exit(0);
			}
		F->trv[i]=msh.kt[e].frent;
		}
	if(F->cs[i]==2)
		{tau=falq_trav_noth(msh,underl[i+1],underl[i]);
		if(tau==-1)
			{fprintf(tmpout,"Unable to find traversed edge\n");
			exit(0);
			}
		F->trv[i]=tau;
		}
	if(F->cs[i]==3)
		{c_inc=mifz_comm_necg(msh,underl[i],underl[i+1]);
		if(c_inc==-1)
			{fprintf(tmpout,"Unable to find common edge\n");
			exit(0);
			}
		F->trv[i]=c_inc;
		}
	if(F->cs[i]==4)
		{tau=falq_trav_noth(msh,underl[i],underl[i+1]);
		if(tau==-1)
			{fprintf(tmpout,"Unable to find traversed edge\n");
			exit(0);
			}
		F->trv[i]=tau;
		}
	}
F->st_grs=nb;
}


void vewk_allo_jovk(int N,float_curve *FC)
{FC->stn =(point *)malloc(N*sizeof(point));
FC->cs       =(int *)malloc(N*sizeof(int));
FC->trv=(int *)malloc(N*sizeof(int));
FC->nd_stat  =(int *)malloc(N*sizeof(int));
FC->nd_idx   =(int *)malloc(N*sizeof(int));
FC->kt_idx   =(int *)malloc(N*sizeof(int));
FC->kt_stat  =(int *)malloc(N*sizeof(int));
}


void lohm_dest_nosr(float_curve *FC)
{free(FC->stn);
free(FC->cs);
free(FC->trv);
free(FC->nd_stat);
free(FC->nd_idx);
free(FC->kt_idx);
free(FC->kt_stat);
}

