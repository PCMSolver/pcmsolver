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
#include "cavity.h"
#include "splinemol.h"
#include "pln_sph.h"
#include "coarsequad.h"



void sejm_allo_bump(int N,span_hz *S)
{int i;
S->nb_children=(int *)malloc(N*sizeof(int));
S->prt     =(int *)malloc(N*sizeof(int));
S->chd   =(int **)malloc(N*sizeof(int*));
for(i=0;i<N;i++)
	S->chd[i]=(int *)malloc(N*sizeof(int));
S->niv_length=(int *)malloc(N*sizeof(int));
S->niv=(int **)malloc(N*sizeof(int*));
for(i=0;i<N;i++)
	S->niv[i]=(int *)malloc(N*sizeof(int));
S->v_value=(int *)malloc(N*sizeof(int));
S->pos=(int *)malloc(N*sizeof(int));
}



void sujt_deal_wozc(int N,span_hz *S)
{int i;
free(S->nb_children);
free(S->prt);
for(i=0;i<N;i++)
	free(S->chd[i]);
free(S->chd);
free(S->niv_length);
for(i=0;i<N;i++)
	free(S->niv[i]);
free(S->niv);
free(S->v_value);
free(S->pos);
}



void nokj_swap_hapk(span_hz *T)
{int nnd,i,j,v,max=40;
span_hz temp;
nnd=T->v_grs;
temp.v_value=(int *)malloc(nnd*sizeof(int));
temp.prt=(int *)malloc(nnd*sizeof(int));
temp.pos=(int *)malloc(nnd*sizeof(int));
temp.nb_children=(int *)malloc(nnd*sizeof(int));
temp.chd=(int **)malloc(nnd*sizeof(int*));
for(i=0;i<nnd;i++)
	temp.chd[i]=(int *)malloc(max*sizeof(int));
for(i=0;i<nnd;i++)
	{v=T->v_value[i];
	temp.prt[v]=T->prt[v];
	temp.pos[v]=T->pos[i];
	temp.nb_children[v]=T->nb_children[v];
	for(j=0;j<temp.nb_children[v];j++)
		temp.chd[v][j]=T->chd[v][j];
	}

for(i=0;i<nnd;i++)
	{T->prt[i]=temp.prt[i];
	T->pos[i]=temp.pos[i];
	T->nb_children[i]=temp.nb_children[i];
	for(j=0;j<temp.nb_children[i];j++)
		T->chd[i][j]=temp.chd[i][j];
	}
free(temp.v_value);
free(temp.prt);
free(temp.pos);
free(temp.nb_children);
for(i=0;i<nnd;i++)
	free(temp.chd[i]);
free(temp.chd);
for(i=0;i<nnd;i++)
	T->v_value[i]=i;
}



double lirn_find_qokc_tecj(prat_main G,int s,span_hz *T)
{int nb,nvt,ned,i,ts1,ts2,idx1,idx2,q,ind;
int n1,n2,e,nch,idx,p,par,L,l,*ED;
double sm,w;
nvt=G.v_grs;
ned=G.k_grs;

nb=1;
T->v_grs=nb;
T->v_value[0]=s;
T->nb_children[s]=0;
T->prt[s]=-1;
T->pos[0]=0;
T->dp=1;
w=0;
ED=(int *)malloc(ned*sizeof(int));
for(i=0;i<ned;i++)
	ED[i]=i;

while(1)
	{sm=9999999.9;
	ind=1;
	for(i=0;i<ned;i++)if(ED[i]!=-1)
		{n1=G.kt[i].str;
		n2=G.kt[i].ter;
		if(G.gew[i]<sm)
			{idx1=-1;	idx2=-1;
			ts1=gonl_arra_govj(T->v_value,T->v_grs,n1,&idx1);
			ts2=gonl_arra_govj(T->v_value,T->v_grs,n2,&idx2);
			if((ts1==1)&&(ts2==0))
				{sm=G.gew[i];
				q=n2;
				e=i;
				idx=idx1;
				ind=2;
				}
			if((ts1==0)&&(ts2==1))
				{sm=G.gew[i];
				q=n1;
				e=i;
				idx=idx2;
				ind=2;
				}
			}
		}
	if(ind==2)
		{par=T->v_value[idx];
		T->v_value[nb]=q;
		T->v_grs=nb+1;
		nch=T->nb_children[par];
		T->chd[par][nch]=q;
		T->prt[q]=par;
		T->nb_children[par]=nch+1;
		T->nb_children[q]=0;
		p=T->pos[idx];
		T->pos[nb]=p+1;
		if(T->pos[nb]+1>T->dp)
			T->dp=T->pos[nb]+1;
		
		w=w+G.gew[e];
		ED[e]=-1;
		nb++;
		}
	if(ind==1)
		break;
	}
free(ED);

for(l=0;l<T->dp;l++)
	{T->niv_length[l]=0;
	for(i=0;i<T->v_grs;i++)
		{q=T->v_value[i];
		p=T->pos[i];
		if(p==l)
			{L=T->niv_length[l];
			T->niv[l][L]=q;
			T->niv_length[l]=L+1;
			}
		}
	}
T->wrz=s;
nokj_swap_hapk(T);
return w;
}



