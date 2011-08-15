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
#include <malloc.h>
#include "cavity.h"
#include "splinemol.h"
#include "pln_sph.h"
#include "coarsequad.h"

 
int wurc_oppo_qijn(manif_tl msh,int n,int s)
{int e[3],res,i,n1,n2;
e[0]=msh.entity[s].frkt;
e[1]=msh.entity[s].sckt;
e[2]=msh.entity[s].trkt;
for(i=0;i<3;i++)
	{n1=msh.kt[e[i]].frvrt;
	n2=msh.kt[e[i]].scvrt;
	if((n1!=n)&&(n2!=n))
		{res=e[i];
		break;
		}
	}
return res;
}



int jedv_chec_tofj(manif_tl msh,int n1,int s,point v_bar)
{int e,na,nb,suc,ps;
pl_eq P;
e=wurc_oppo_qijn(msh,n1,s);
na=msh.kt[e].frvrt;
nb=msh.kt[e].scvrt;
vowg_plan_vowt(msh.knot[na],msh.knot[nb],msh.knot[n1],&P);
ps=balt_posi_wurp(msh.knot[n1],v_bar,P);
if(ps==1)
	suc=SUCCESS;
else
	suc=FAILURE;
return suc;
}



int caqt_excl_nukh(manif_tl msh,int n1,int n2,int *el,int max_exc)
{int nb,*temp1,*temp2,val1,val2,i,k,ts,dummy;
temp1=(int *)malloc(2*msh.increl[n1].val*sizeof(int));
temp2=(int *)malloc(2*msh.increl[n2].val*sizeof(int));
val1=cepj_inci_jeqt(msh,n1,temp1);
val2=cepj_inci_jeqt(msh,n2,temp2);
k=0;
for(i=0;i<val1;i++)
	{ts=gonl_arra_govj(temp2,val2,temp1[i],&dummy);
	if(ts==0)
		{if(k>=max_exc)
			{fprintf(tmpout,"max_exc is attained\n");
			exit(0);
			}
		el[k]=temp1[i];
		k++;
		}
	}
free(temp1);
free(temp2);
nb=k;
return nb;
}



int lufj_cons_retw(manif_tl msh,int n1,int n2,point v_bar)
{int suc=SUCCESS,val,*el,i,ts,s,max_exc=100;
el=(int *)malloc(max_exc*sizeof(int));

val=caqt_excl_nukh(msh,n1,n2,el,max_exc);
for(i=0;i<val;i++)
	{s=el[i];
	ts=jedv_chec_tofj(msh,n1,s,v_bar);
	if(ts==FAILURE)
		{suc=FAILURE;
		break;
		}
	}

val=caqt_excl_nukh(msh,n2,n1,el,max_exc);
for(i=0;i<val;i++)
	{s=el[i];
	ts=jedv_chec_tofj(msh,n2,s,v_bar);
	if(ts==FAILURE)
		{suc=FAILURE;
		break;
		}
	}
free(el);
return suc;
}



int local_node_integrity_vil(manif_tl msh,int s)
{int ts=1,val,i,f,n1,n2;
val=msh.increl[s].val;
for(i=0;i<val;i++)
	{f=msh.increl[s].inc[i];
	n1=msh.kt[f].frvrt;
	n2=msh.kt[f].scvrt;
	if((n1!=s)&&(n2!=s))
		{ts=0;
		break;
		}
	}
if(ts==0)
	{fprintf(tmpout,"node=%d with bad integrity\n",s);
	exit(0);
	}
return ts;
}


void picv_chec_qedm(manif_tl msh,int p,int q)
{int n1,n2,n3,n_loc[3],ts1,ts2,ts3,dummy;
n1=msh.entity[p].frvrt;
n2=msh.entity[p].scvrt;
n3=msh.entity[p].thvrt;
n_loc[0]=msh.entity[q].frvrt;
n_loc[1]=msh.entity[q].scvrt;
n_loc[2]=msh.entity[q].thvrt;
ts1=gonl_arra_govj(n_loc,3,n1,&dummy);
ts2=gonl_arra_govj(n_loc,3,n2,&dummy);
ts3=gonl_arra_govj(n_loc,3,n3,&dummy);
if((ts1==1)&&(ts2==1)&&(ts3==1))
	{fprintf(tmpout,"element permutation\n");
	exit(0);
	}
}


void gewn_chec_qamk(manif_tl msh)
{int nnd,nel,ned,s;

nnd=msh.n_grs;
nel=msh.e_grs;
ned=msh.k_grs;


for(s=0;s<nel;s++)
	sokm_loca_vazn(msh,s);

for(s=0;s<ned;s++)
	cugt_loca_pogt(msh,s);

for(s=0;s<nnd;s++)
	local_node_integrity_vil(msh,s);

}

 
void mokq_chec_nukb(manif_tl msh)
{int nnd,nel,ned,s,cld;

nnd=msh.n_grs;
nel=msh.e_grs;
ned=msh.k_grs;


for(s=0;s<nel;s++)
	sokm_loca_vazn(msh,s);

cld=0;
for(s=0;s<ned;s++)
	{darj_loca_wejf(msh,s);
	if(msh.kt[s].scent==-1)
		cld=1;
	}
if(cld==0)	fprintf(tmpout,"manifold is closed\n");
if(cld==1)	fprintf(tmpout,"manifold is open\n");
fprintf(tmpout,"Good manifold integrity without incre\n");
}



void rapv_coll_moks(manif_tl *msh,int *el,int *e,
point P,int *map_node,int *map_edge)
{int n1,n2,*temp,ned,*map1,*map2,new_e0,s,i,m;
ned=msh->k_grs;
map1=(int *)malloc(ned*sizeof(int));
dopr_fuse_tenv(msh,e[1],e[2],el[3],el[4],e[3],e[4],el[5],el[6],map1);
new_e0=map1[e[0]];

if(msh->kt[new_e0].frvrt<msh->kt[new_e0].scvrt)
	{n1=msh->kt[new_e0].frvrt;
	n2=msh->kt[new_e0].scvrt;
	}
else
	{n2=msh->kt[new_e0].frvrt;
	n1=msh->kt[new_e0].scvrt;
	}

roph_coll_fact(msh,n1,n2,P);
fuwr_repl_wejf(msh,n1,n2);
zudp_disc_duwc(msh,n2,map_node);

temp=(int *)malloc(2*sizeof(int));
temp[0]=el[1];
temp[1]=el[2];
vegs_disc_lepq(msh,temp,2);
free(temp);

map2=(int *)malloc(ned*sizeof(int));
m=msh->k_grs;
topr_disc_nufr(msh,new_e0,map2);
for(i=0;i<ned;i++)
	{s=map1[i];
	if(s==-1)
		map_edge[i]=-1;
	else
		map_edge[i]=map2[s];
	}
free(map1);
free(map2);
}



void newq_rema_furw(manif_tl msh,int e,int E,int *f,int *F)
{int i,k,ed[3],r,s;
ed[0]=msh.entity[E].frkt;
ed[1]=msh.entity[E].sckt;
ed[2]=msh.entity[E].trkt;
k=0;
for(i=0;i<3;i++)
	{r=msh.kt[ed[i]].frent;
	s=msh.kt[ed[i]].scent;
	if(ed[i]!=e)
		{f[k]=ed[i];
		if(r!=E)	F[k]=r;
		else		F[k]=s;
		k++;
		if(k==2)
			break;
		}
	}
}



void pajc_coll_qetv(manif_tl *msh,int f,point P,
int *map_node,int *map_edge)
{int *e,*el,*g,*G;
e=(int *)malloc(5*sizeof(int));
el=(int *)malloc(7*sizeof(int));
e[0]=f;
el[1]=msh->kt[f].frent;
el[2]=msh->kt[f].scent;

g=(int *)malloc(2*sizeof(int));
G=(int *)malloc(2*sizeof(int));
newq_rema_furw(*msh,f,el[1],g,G);
if(g[0]<g[1])
	{e[1]=g[0];	el[3]=G[0];
	e[2]=g[1];	el[4]=G[1];
	}
else
	{e[1]=g[1];	el[3]=G[1];
	e[2]=g[0];	el[4]=G[0];
	}
newq_rema_furw(*msh,f,el[2],g,G);
if(g[0]<g[1])
	{e[3]=g[0];	el[5]=G[0];
	e[4]=g[1];	el[6]=G[1];
	}
else
	{e[4]=g[0];	el[6]=G[0];
	e[3]=g[1];	el[5]=G[1];
	}
free(g);
free(G);
rapv_coll_moks(msh,el,e,P,map_node,map_edge);
free(e);
free(el);
}




