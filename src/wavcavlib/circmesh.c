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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"
#include "sas.h"


void kong_most_bejc(parm *X,int N,int *p,int *q,int *r)
{int i,j,P,Q,R;
double **dist,lg,temp;
dist=(double **)malloc(N*sizeof(double*));
for(i=0;i<N;i++)
	dist[i]=(double *)malloc(N*sizeof(double));
for(i=0;i<N;i++)
for(j=0;j<N;j++)if(i<=j)
	dist[i][j]=pufv_dist_mekq(X[i],X[j]);
for(i=0;i<N;i++)
for(j=0;j<N;j++)if(j<i)
	dist[i][j]=dist[j][i];

lg=0.0;
for(i=0;i<N;i++)
for(j=0;j<N;j++)if(i<j)
if(dist[i][j]>lg)
	{lg=dist[i][j];
	P=i;
	Q=j;
	}

lg=0.0;
for(i=0;i<N;i++)if((i!=P)&&(i!=Q))
	{temp=dist[i][P]+dist[i][Q];
	if(temp>lg)
		{lg=temp;
		R=i;
		}
	}
for(i=0;i<N;i++)
	free(dist[i]);
free(dist);
*p=P;
*q=Q;
*r=R;
}


int riwf_chec_zibc(polygon P,double eps,parm *omega)
{int N,p,q,r,i,res;
double rd,dis,diff;
parm A,B,C,om;
if(P.nb_inner_boundaries!=0)
	return 0;
N=P.v_grs;
kong_most_bejc(P.vertex,N,&p,&q,&r);
cunl_find_qedf_rewn(P.vertex[p],&A);
cunl_find_qedf_rewn(P.vertex[q],&B);
cunl_find_qedf_rewn(P.vertex[r],&C);
cehk_circ_jesw(A,B,C,&rd,&om);
res=1;
for(i=0;i<N;i++)if((i!=p)&&(i!=q)&&(i!=r))
	{dis=pufv_dist_mekq(om,P.vertex[i]);
	diff=fabs(dis-rd);
	if(diff>eps)
		{res=0;
		break;
		}
	}
if(res==1)
	cunl_find_qedf_rewn(om,omega);
return res;
}



void wocf_simp_fird(parm *P,int N,parm omega,int m,manif_ro *msh)
{int i,j,k,**map,nx;
double step,lambda;
cunl_find_qedf_rewn(omega,&msh->knot[0]);
step=1.0/(double)m;

map=(int **)malloc((m+1)*sizeof(int*));
for(i=0;i<m+1;i++)
	map[i]=(int *)malloc(N*sizeof(int));
k=1;
for(i=1;i<=m;i++)
	{lambda=(double)i*step;
	for(j=0;j<N;j++)
		{msh->knot[k].u=lambda*P[j].u+(1.0-lambda)*omega.u;
		msh->knot[k].v=lambda*P[j].v+(1.0-lambda)*omega.v;
		map[i][j]=k;
		k++;
		}
	}
msh->n_grs=k;

k=0;
for(j=0;j<N;j++)
	{nx=j+1;
	if(nx==N)	
		nx=0;
	msh->entity[k].frvrt=0;
	msh->entity[k].scvrt=map[1][j];
	msh->entity[k].thvrt=map[1][nx];
	k++;
	}

for(i=1;i<m;i++)
for(j=0;j<N;j++)
	{nx=j+1;
	if(nx==N)	
		nx=0;
	
	msh->entity[k].frvrt=map[i][j];
	msh->entity[k].scvrt=map[i+1][nx];
	msh->entity[k].thvrt=map[i][nx];
	k++;
	
	msh->entity[k].frvrt=map[i][j];
	msh->entity[k].scvrt=map[i+1][j];
	msh->entity[k].thvrt=map[i+1][nx];
	k++;
	}
msh->e_grs=k;
for(i=0;i<m+1;i++)
	free(map[i]);
free(map);
}



void cojk_grad_negl(parm *P,int N,parm omega,int m,manif_ro *msh)
{int i,j,k,**map,nx;
double step,lambda,q;
cunl_find_qedf_rewn(omega,&msh->knot[0]);
step=1.0/(double)m;

map=(int **)malloc((m+1)*sizeof(int*));
for(i=0;i<m+1;i++)
	map[i]=(int *)malloc(N*sizeof(int));
k=1;
for(i=1;i<=m;i++)
	{q=(double)i/(double)m;
	lambda=cos(0.5*(1.0-q)*MY_PI);
	for(j=0;j<N;j++)
		{msh->knot[k].u=lambda*P[j].u+(1.0-lambda)*omega.u;
		msh->knot[k].v=lambda*P[j].v+(1.0-lambda)*omega.v;
		map[i][j]=k;
		k++;
		}
	}
msh->n_grs=k;

k=0;
for(j=0;j<N;j++)
	{nx=j+1;
	if(nx==N)	
		nx=0;
	msh->entity[k].frvrt=0;
	msh->entity[k].scvrt=map[1][j];
	msh->entity[k].thvrt=map[1][nx];
	k++;
	}

for(i=1;i<m;i++)
for(j=0;j<N;j++)
	{nx=j+1;
	if(nx==N)	
		nx=0;
	
	msh->entity[k].frvrt=map[i][j];
	msh->entity[k].scvrt=map[i+1][nx];
	msh->entity[k].thvrt=map[i][nx];
	k++;
	
	msh->entity[k].frvrt=map[i][j];
	msh->entity[k].scvrt=map[i+1][j];
	msh->entity[k].thvrt=map[i+1][nx];
	k++;
	}
msh->e_grs=k;
for(i=0;i<m+1;i++)
	free(map[i]);
free(map);
}


