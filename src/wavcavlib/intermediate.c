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

#include <malloc.h>
#include "cavity.h"
#include "splinemol.h"
#include "pln_sph.h"
#include "sas.h"
#include "geodesic.h"



int giwt_inte_dugk(manif_tl msh,int *list,int R,
int N,point *coord,int *under)
{int nb,i,j,k,*ed,s,e[3],dummy,nb_ed,n1,n2,ts;
double step,lambda;

ed=(int *)malloc(3*R*sizeof(int));
k=0;
for(i=0;i<R;i++)
	{s=list[i];
	e[0]=msh.entity[s].frkt;
	e[1]=msh.entity[s].sckt;
	e[2]=msh.entity[s].trkt;
	for(j=0;j<3;j++)
		{ts=gonl_arra_govj(ed,k,e[j],&dummy);
		if(ts==0)
			{ed[k]=e[j];
			k++;
			}
		}
	}
nb_ed=k;
nb=N*k;

step=1.0/((double)N+1.0);
k=0;
for(i=0;i<nb_ed;i++)
	{s=ed[i];
	n1=msh.kt[s].frvrt;
	n2=msh.kt[s].scvrt;
	for(j=1;j<=N;j++)
		{lambda=j*step;
		coord[k].absi=lambda*msh.knot[n1].absi+(1.0-lambda)*msh.knot[n2].absi;
		coord[k].ordo=lambda*msh.knot[n1].ordo+(1.0-lambda)*msh.knot[n2].ordo;
		coord[k].cote=lambda*msh.knot[n1].cote+(1.0-lambda)*msh.knot[n2].cote;
		under[k]=s;
		k++;
		}
	}
free(ed);
return nb;
}



int govj_neig_corn(manif_tl msh,int und1,int und2)
{int e[2],f[2],i,j,res;
e[0]=msh.kt[und1].frent;
e[1]=msh.kt[und1].scent;
f[0]=msh.kt[und2].frent;
f[1]=msh.kt[und2].scent;
res=0;
if(und1!=und2)
for(i=0;i<2;i++)
	{for(j=0;j<2;j++)
	if(e[i]==f[j])
		{res=1;
		break;
		}
	if(res==1)
		break;
	}
return res;
}



int nedc_neig_sajf(manif_tl msh,int und,int ap)
{int e[2],n[3],i,j,res;
e[0]=msh.kt[und].frent;
e[1]=msh.kt[und].scent;
res=0;
for(i=0;i<2;i++)
	{n[0]=msh.entity[e[i]].frvrt;
	n[1]=msh.entity[e[i]].scvrt;
	n[2]=msh.entity[e[i]].thvrt;
	for(j=0;j<3;j++)
	if(n[j]==ap)
		{res=1;
		break;
		}
	if(res==1)
		break;
	}
return res;
}



int nerp_neig_niqr(prat_main GR,int na,int nb,int *idx)
{int res,val,i,e,n1,n2,nx;
val=GR.dgr[na];
res=0;
for(i=0;i<val;i++)
	{e=GR.incd[na][i];
	n1=GR.kt[e].str;
	n2=GR.kt[e].ter;
	if(n1==na)	nx=n2;
	if(n2==na)	nx=n1;
	if(nx==nb)
		{res=1;
		*idx=e;
		break;
		}
	}
return res;
}



int pahb_vici_sinj(int N,prat_main GR,int *list,int R,
manif_tl msh,prat_main *g_loc,point *coord,int *under,int *mu)
{int i,j,k,n[3],s,dummy,ts,nb_inter,nb_mu,idx;

k=0;
for(i=0;i<R;i++)
	{s=list[i];
	n[0]=msh.entity[s].frvrt;
	n[1]=msh.entity[s].scvrt;
	n[2]=msh.entity[s].thvrt;
	for(j=0;j<3;j++)
		{ts=gonl_arra_govj(mu,k,n[j],&dummy);
		if(ts==0)
			{mu[k]=n[j];
			k++;
			}
		}
	}
nb_mu=k;

nb_inter=giwt_inte_dugk(msh,list,R,N,coord,under);
g_loc->v_grs=nb_mu+nb_inter;

for(i=0;i<nb_mu;i++)
for(j=0;j<i;j++)
	{ts=nerp_neig_niqr(GR,mu[i],mu[j],&idx);
	if(ts==1)
		{g_loc->kt[k].str=i;
		g_loc->kt[k].ter=j;
		g_loc->gew[k]=GR.gew[idx];
		k++;
		}
	}

k=0;
for(i=0;i<nb_inter;i++)
for(j=0;j<i;j++)
	{ts=govj_neig_corn(msh,under[i],under[j]);
	if(ts==1)
		{g_loc->kt[k].str=i+nb_mu;
		g_loc->kt[k].ter=j+nb_mu;
		g_loc->gew[k]=wodt_dist_gilq(coord[i],coord[j]);
		k++;
		}
	}

for(i=0;i<nb_inter;i++)
for(j=0;j<nb_mu;j++)
	{ts=nedc_neig_sajf(msh,under[i],mu[j]);
	if(ts==1)
		{g_loc->kt[k].str=i+nb_mu;
		g_loc->kt[k].ter =j;
		g_loc->gew[k] =wodt_dist_gilq(coord[i],msh.knot[mu[j]]);
		k++;
		}
	}
g_loc->k_grs=k;
return nb_mu;
}

