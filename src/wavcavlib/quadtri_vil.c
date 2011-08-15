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
#include <malloc.h>
#include <stdlib.h>
#include "cavity.h"
#include "splinemol.h" 
#include "pln_sph.h"
#include "coarsequad.h"


void gapd_remo_purw(span_hz *T,int s,int z)
{int i,q,L;
L=T->niv_length[s];
for(i=0;i<L;i++)
if(T->niv[s][i]==z)
	{q=i;
	break;
	}
for(i=q;i<L-1;i++)
	T->niv[s][i]=T->niv[s][i+1];
T->niv_length[s]=L-1;
}


void qetk_trim_zogf(span_hz *T,int z)
{int nvt,l,d,L,i,p,q,nb;
if(z==T->wrz)
	{T->v_grs=0;
	T->dp=0;
	}
else
	{nvt=T->v_grs;
	l=T->pos[z];
	gapd_remo_purw(T,l,z);
	
	p=T->prt[z];
	nb=T->nb_children[p];
	for(i=0;i<nb;i++)
	if(T->chd[p][i]==z)
		{q=i;
		break;
		}
	for(i=q;i<nb-1;i++)
		T->chd[p][i]=T->chd[p][i+1];
	T->nb_children[p]=nb-1;
	
	d=T->dp;
	L=T->niv_length[d-1];
	if(L==0)
		T->dp=d-1;
	
	T->v_grs=nvt-1;
	}
}


int dorv_comp_qomr(manif_tl msh,int s,int r)
{int res,*nd,ts,dummy;
nd=(int *)malloc(3*sizeof(int));
nd[0]=msh.entity[r].frvrt;
nd[1]=msh.entity[r].scvrt;
nd[2]=msh.entity[r].thvrt;
ts=gonl_arra_govj(nd,3,msh.entity[s].frvrt,&dummy);
if(ts==0)	res=msh.entity[s].frvrt;
ts=gonl_arra_govj(nd,3,msh.entity[s].scvrt,&dummy);
if(ts==0)	res=msh.entity[s].scvrt;
ts=gonl_arra_govj(nd,3,msh.entity[s].thvrt,&dummy);
if(ts==0)	res=msh.entity[s].thvrt;
free(nd);
return res;
}


void kepg_mela_wekj(manif_tl msh,fajor_sion3D *QUAD,
int u,int w,point omega,int max_nnd,int max_nel)
{int h,k,m,nx;
h=QUAD->e_grs;
if(h>=max_nel)
	{fprintf(tmpout,"max_nel=%d is obtained\n",max_nel);
	exit(0);
	}
k=QUAD->n_grs;
m=dorv_comp_qomr(msh,u,w);
if(k>=max_nnd)
	{fprintf(tmpout,"1-max_nnd=%d is obtained\n",max_nnd);
	exit(0);
	}
QUAD->knot[k].absi=msh.knot[m].absi;
QUAD->knot[k].ordo=msh.knot[m].ordo;
QUAD->knot[k].cote=msh.knot[m].cote;
QUAD->elem[h].frvrt=k;
k++;
nx=kuts_next_kolv(msh.entity[u],m);
m=nx;
if(k>=max_nnd)
	{fprintf(tmpout,"2-max_nnd=%d is obtained\n",max_nnd);
	exit(0);
	}
QUAD->knot[k].absi=msh.knot[m].absi;
QUAD->knot[k].ordo=msh.knot[m].ordo;
QUAD->knot[k].cote=msh.knot[m].cote;
QUAD->elem[h].scvrt=k;
k++;
if(k>=max_nnd)
	{fprintf(tmpout,"3-max_nnd=%d is obtained\n",max_nnd);
	exit(0);
	}
QUAD->knot[k].absi=omega.absi;
QUAD->knot[k].ordo=omega.ordo;
QUAD->knot[k].cote=omega.cote;
QUAD->elem[h].thvrt=k;
k++;
nx=kuts_next_kolv(msh.entity[u],m);
m=nx;
if(k>=max_nnd)
	{fprintf(tmpout,"4-max_nnd=%d is obtained\n",max_nnd);
	exit(0);
	}
QUAD->knot[k].absi=msh.knot[m].absi;
QUAD->knot[k].ordo=msh.knot[m].ordo;
QUAD->knot[k].cote=msh.knot[m].cote;
QUAD->elem[h].ftvrt=k;
k++;
QUAD->n_grs=k;
QUAD->e_grs=h+1;
}


int rilk_comm_ditl(manif_tl msh,int s,int t,int *res)
{int nb,n[3],m[3],i,j,k;
n[0]=msh.entity[s].frvrt;
n[1]=msh.entity[s].scvrt;
n[2]=msh.entity[s].thvrt;
m[0]=msh.entity[t].frvrt;
m[1]=msh.entity[t].scvrt;
m[2]=msh.entity[t].thvrt;
k=0;
for(i=0;i<3;i++)
for(j=0;j<3;j++)
if(n[i]==m[j])
	{res[k]=n[i];
	k++;
	}
nb=k;
return nb;
}


void decs_quad_sibl(span_hz T,int u,manif_tl *msh,
fajor_sion3D *QUAD,int max_nnd,int max_nel)
{int v,w,c,nnd;
int p,n1,n2,n3;
point omega;

p=T.prt[u];
if(T.chd[p][0]==u)	v=T.chd[p][1];
if(T.chd[p][1]==u)	v=T.chd[p][0];

w=T.prt[u];
n1=msh->entity[w].frvrt;
n2=msh->entity[w].scvrt;
n3=msh->entity[w].thvrt;
omega.absi=(msh->knot[n1].absi+msh->knot[n2].absi+msh->knot[n3].absi)/3.0;
omega.ordo=(msh->knot[n1].ordo+msh->knot[n2].ordo+msh->knot[n3].ordo)/3.0;
omega.cote=(msh->knot[n1].cote+msh->knot[n2].cote+msh->knot[n3].cote)/3.0;

kepg_mela_wekj(*msh,QUAD,u,w,omega,max_nnd,max_nel);
kepg_mela_wekj(*msh,QUAD,v,w,omega,max_nnd,max_nel);
rilk_comm_ditl(*msh,u,v,&c);
nnd=msh->n_grs;
if(nnd>=max_nnd)
	{fprintf(tmpout,"5-max_nnd=%d is obtained\n",max_nnd);
	exit(0);
	}
msh->knot[nnd].absi=omega.absi;
msh->knot[nnd].ordo=omega.ordo;
msh->knot[nnd].cote=omega.cote;
if(msh->entity[w].frvrt==c)	msh->entity[w].frvrt=nnd;
if(msh->entity[w].scvrt==c)	msh->entity[w].scvrt=nnd;
if(msh->entity[w].thvrt==c)	msh->entity[w].thvrt=nnd;
msh->n_grs=nnd+1;
}


void lomr_quad_gudw(span_hz T,int u,manif_tl msh,
fajor_sion3D *QUAD,int max_nnd,int max_nel)
{int w,c,h,k,m,nx;
w=T.prt[u];
h=QUAD->e_grs;
k=QUAD->n_grs;
if(h>=max_nel)
	{fprintf(tmpout,"max_nel=%d is attained\n",max_nel);
	exit(0);
	}

c=dorv_comp_qomr(msh,u,w);
if(k>=max_nnd)
	{fprintf(tmpout,"1-max_nnd is attained (%d)\n",max_nnd);
	exit(0);
	}
QUAD->knot[k].absi=msh.knot[c].absi;
QUAD->knot[k].ordo=msh.knot[c].ordo;
QUAD->knot[k].cote=msh.knot[c].cote;
QUAD->elem[h].frvrt=k;
k++;

nx=kuts_next_kolv(msh.entity[u],c);
m=nx;
if(k>=max_nnd)
	{fprintf(tmpout,"2-max_nnd is attained (%d)\n",max_nnd);
	exit(0);
	}
QUAD->knot[k].absi=msh.knot[m].absi;
QUAD->knot[k].ordo=msh.knot[m].ordo;
QUAD->knot[k].cote=msh.knot[m].cote;
QUAD->elem[h].scvrt=k;
k++;

nx=kuts_next_kolv(msh.entity[w],m);
m=nx;
if(k>=max_nnd)
	{fprintf(tmpout,"3-max_nnd is attained (%d)\n",max_nnd);
	exit(0);
	}
QUAD->knot[k].absi=msh.knot[m].absi;
QUAD->knot[k].ordo=msh.knot[m].ordo;
QUAD->knot[k].cote=msh.knot[m].cote;
QUAD->elem[h].thvrt=k;
k++;

nx=kuts_next_kolv(msh.entity[w],m);
m=nx;
if(k>=max_nnd)
 	{fprintf(tmpout,"4-max_nnd is attained (%d)\n",max_nnd);
	exit(0);
	}
QUAD->knot[k].absi=msh.knot[m].absi;
QUAD->knot[k].ordo=msh.knot[m].ordo;
QUAD->knot[k].cote=msh.knot[m].cote;
QUAD->elem[h].ftvrt=k;
k++;
QUAD->n_grs=k;
QUAD->e_grs=h+1;
}


void vocm_disp_qosh(span_hz S)
{int i,j,v;
fprintf(tmpout,"Number of vertices=%d\n",S.v_grs);
for(i=0;i<S.v_grs;i++)
	{v=S.v_value[i];
	fprintf(tmpout,"vertex=%d   parent=%d  children=",v,S.prt[v]);
	for(j=0;j<S.nb_children[v];j++)
		fprintf(tmpout,"%d  ",S.chd[v][j]);
	fprintf(tmpout,"\n");
	}
fprintf(tmpout,"Depth=%d\n",S.dp);
for(i=0;i<S.dp;i++)
	{fprintf(tmpout,"Level=%d  members=",i);
	for(j=0;j<S.niv_length[i];j++)
		fprintf(tmpout,"%d  ",S.niv[i][j]);
	fprintf(tmpout,"\n");
	}
fprintf(tmpout,"----------------------------\n");
}


void vogn_allo_cusr(int nnd,int ned,int max_inc,prat_main *G)
{int i;
G->kt=(bingo *)malloc(ned*sizeof(bingo));
G->dgr=(int *)malloc(nnd*sizeof(int));
G->gew=(double *)malloc(ned*sizeof(double));
G->incd=(int **)malloc(nnd*sizeof(int*));
for(i=0;i<nnd;i++)
	G->incd[i]=(int *)malloc(max_inc*sizeof(int));
}


void conw_dest_vojk(int nnd,prat_main *G)
{int i;
free(G->kt);
free(G->dgr);
free(G->gew);
for(i=0;i<nnd;i++)
	free(G->incd[i]);
free(G->incd);
}


void maqb_comp_sogf(manif_tl msh,manif_tl *out)
{int nnd,nel,ned,i; 
nnd=msh.n_grs;
nel=msh.e_grs;
ned=msh.k_grs;
for(i=0;i<nnd;i++)
	{out->knot[i].absi=msh.knot[i].absi;
	out->knot[i].ordo=msh.knot[i].ordo;
	out->knot[i].cote=msh.knot[i].cote;
	}
for(i=0;i<nel;i++)
	{out->entity[i].frvrt=msh.entity[i].frvrt;
	out->entity[i].scvrt=msh.entity[i].scvrt;
	out->entity[i].thvrt=msh.entity[i].thvrt;
	out->entity[i].frkt=msh.entity[i].frkt;
	out->entity[i].sckt=msh.entity[i].sckt;
	out->entity[i].trkt=msh.entity[i].trkt;
	}
for(i=0;i<ned;i++)
	{out->kt[i].frvrt=msh.kt[i].frvrt;
	out->kt[i].scvrt=msh.kt[i].scvrt;
	out->kt[i].frent=msh.kt[i].frent;
	out->kt[i].scent=msh.kt[i].scent;
	}
out->n_grs=nnd;
out->e_grs=nel;
out->k_grs=ned;
}


int dosc_coor_licf(point *X,int n,
point Y,double eps,int *id)
{int res,i,ts;
res=0;
for(i=0;i<n;i++)
	{ts=gect_tole_husn(Y,X[i],eps);
	if(ts==1)
		{res=1;
		*id=i;
		break;
		}
	}
return res;
}


