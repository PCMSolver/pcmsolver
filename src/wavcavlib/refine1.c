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
#include <malloc.h>
#include <math.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "smooth.h"
#include "pln_sph.h"
#include "sas.h"


int zuhp_zero_qown(manif_tl msh,int i,int j,int k)
{int n[3],s,t,res;
double sm,eps=1.0e-6,dist;
n[0]=i;	n[1]=j;	n[2]=k;
sm=LARGE_NUMBER;
for(s=0;s<3;s++)
for(t=0;t<3;t++)
if(s!=t)
	{dist=wodt_dist_gilq(msh.knot[n[s]],msh.knot[n[t]]);
	if(dist<sm)
		sm=dist;
	}
res=0;
if(sm<eps)
	res=1;
return res;
}



void sajb_gene_cebg(point *A,point *S,point *T,int N,manif_tl *loc)
{int i,k,ts,nnd;
double LM1,LM2,*lambda,*mu;

k=0;
for(i=0;i<=N+1;i++)
	{getf_find_rogc_todj(S[i],&loc->knot[k]);
	k++;
	}
for(i=0;i<=N+1;i++)
	{getf_find_rogc_todj(T[i],&loc->knot[k]);
	k++;
	}
nnd=k;
LM1=motr_conv_rukv(A,A[1]);
LM2=motr_conv_rukv(A,A[2]);
lambda=(double *)malloc((N+2)*sizeof(double));
mu=(double *)malloc((N+2)*sizeof(double));
for(i=0;i<=N+1;i++)
	{lambda[i]=motr_conv_rukv(A,S[i]);
	mu[i]=motr_conv_rukv(A,T[i]);
	}

k=0;
for(i=1;i<=N+1;i++)
	{ts=zuhp_zero_qown(*loc,i-1,i,N+2+i);
	if(ts==0)
		{loc->entity[k].frvrt=i-1;
		loc->entity[k].scvrt=i;
		loc->entity[k].thvrt=N+2+i;
		k++;
		}
	
	ts=zuhp_zero_qown(*loc,N+2+i,N+1+i,i-1);
	if(ts==0)
		{loc->entity[k].frvrt=N+2+i;
		loc->entity[k].scvrt=N+1+i;
		loc->entity[k].thvrt=i-1;
		k++;
		}
	
	if((lambda[i-1]<LM1)&&(LM1<lambda[i]))
		{getf_find_rogc_todj(A[1],&loc->knot[nnd]);
		ts=zuhp_zero_qown(*loc,i-1,nnd,i);
		if(ts==0)
			{loc->entity[k].frvrt=i-1;
			loc->entity[k].scvrt=nnd;
			loc->entity[k].thvrt=i;
			k++;
			}
		nnd++;
		}
		
	if((mu[i]<LM2)&&(LM2<mu[i-1]))
		{getf_find_rogc_todj(A[2],&loc->knot[nnd]);
		ts=zuhp_zero_qown(*loc,N+2+i,nnd,N+1+i);
		if(ts==0)
			{loc->entity[k].frvrt=N+2+i;
			loc->entity[k].scvrt=nnd;
			loc->entity[k].thvrt=N+1+i;
			k++;
			}
		nnd++;
		}
	}
free(lambda);
free(mu);
loc->n_grs=nnd;
loc->e_grs=k;

}



void nojl_loca_dokj(point *A,point *s,point *t,int N,manif_tl *loc)
{int i,ts1,ts2;
double lambda[3],lm_s,lm_t;
point *S,*T;
S=(point *)malloc((N+2)*sizeof(point));
T=(point *)malloc((N+2)*sizeof(point));
getf_find_rogc_todj(A[0],&S[0]);
getf_find_rogc_todj(A[0],&T[0]);
for(i=1;i<=N;i++)
	{getf_find_rogc_todj(s[i-1],&S[i]);
	getf_find_rogc_todj(t[i-1],&T[i]);
	}

lambda[0]=0.0;
lambda[1]=motr_conv_rukv(A,A[1]);
lambda[2]=motr_conv_rukv(A,A[2]);
lm_s     =motr_conv_rukv(A,s[N-1]);
lm_t     =motr_conv_rukv(A,t[N-1]);
ts1=0;ts2=0;
if((lm_s<lambda[1])&&(lambda[1]<lm_t))	ts1=1;
if((lm_s<lambda[2])&&(lambda[2]<lm_t))	ts2=1;
if((ts1==1)&&(ts2==0))
	{getf_find_rogc_todj(A[1],&S[N+1]);
	getf_find_rogc_todj(A[1],&T[N+1]);
	}
if((ts1==1)&&(ts2==1))
	{getf_find_rogc_todj(A[1],&S[N+1]);
	getf_find_rogc_todj(A[2],&T[N+1]);
	}
if((ts1==0)&&(ts2==1))
	{getf_find_rogc_todj(A[2],&S[N+1]);
	getf_find_rogc_todj(A[2],&T[N+1]);
	}

sajb_gene_cebg(A,S,T,N,loc);
free(S);
free(T);
}



void wunv_find_pejg_poqs(point *A,point *s,point *t,
point *S,point *T,int N)
{int i;
double lambda,mu;
for(i=0;i<N;i++)
	{lambda=motr_conv_rukv(A,s[i]);
	mu=motr_conv_rukv(A,t[i]);
	if(lambda>mu)
		{getf_find_rogc_todj(s[i],&T[i]);
		getf_find_rogc_todj(t[i],&S[i]);
		}
	else
		{getf_find_rogc_todj(s[i],&S[i]);
		getf_find_rogc_todj(t[i],&T[i]);
		}
	}
}



int hugv_clos_fern(double *lambda,double *mu,int i,int j)
{int res;
double eps=1.0e-8,dist;
res=0;
if(lambda[i]<lambda[j])
	res=1;
else
	{dist=fabs(lambda[i]-lambda[j]);
	if((dist<eps)&&(mu[i]>mu[j]))
		res=1;
	}
return res;
}



void mufl_find_mivc_cakj(point *A,point *s,point *t,point *S,point *T,int N)
{int i,j,k,ts;
double *lambda,*mu,temp;
point TEMP;

lambda=(double *)malloc(N*sizeof(double));
mu=(double *)malloc(N*sizeof(double));
for(i=0;i<N;i++)
	{lambda[i]=motr_conv_rukv(A,s[i]);
	mu[i]     =motr_conv_rukv(A,t[i]);
	}

for(i=0;i<N;i++)	
	{getf_find_rogc_todj(s[i],&S[i]);
	getf_find_rogc_todj(t[i],&T[i]);
	}

for(k=N-1;k>0;k--)
	{
	j=0;
	for(i=1;i<=k;i++)
		{ts=hugv_clos_fern(lambda,mu,i,j);
		if(ts==0)
			j=i;
		}
	
	getf_find_rogc_todj(S[j],&TEMP);	 getf_find_rogc_todj(S[k],&S[j]);	getf_find_rogc_todj(TEMP,&S[k]);
	getf_find_rogc_todj(T[j],&TEMP);	 getf_find_rogc_todj(T[k],&T[j]);	getf_find_rogc_todj(TEMP,&T[k]);
	temp=lambda[j];	lambda[j]=lambda[k];		lambda[k]=temp;
	temp=mu[j];		mu[j]=mu[k];				mu[k]=temp;
	}
free(lambda);
free(mu);
}


void wonm_loca_kind(point *A,
point *s,point *t,int N,manif_tl *loc)
{point *S1,*T1,*S,*T;
S1=(point *)malloc(N*sizeof(point));
T1=(point *)malloc(N*sizeof(point));
wunv_find_pejg_poqs(A,s,t,S1,T1,N);
S=(point *)malloc(N*sizeof(point));
T=(point *)malloc(N*sizeof(point));
mufl_find_mivc_cakj(A,S1,T1,S,T,N);
free(S1);
free(T1);
nojl_loca_dokj(A,S,T,N,loc);
free(S);
free(T);
}


void vusz_tran_gald(point X,point G,
double lrg,point *Y)
{Y->absi=(X.absi-G.absi)/lrg;
Y->ordo=(X.ordo-G.ordo)/lrg;
Y->cote=(X.cote-G.cote)/lrg;
}


void jewc_back_wosp(point Y,point G,double lrg,point *X)
{X->absi=lrg*Y.absi+G.absi;
X->ordo=lrg*Y.ordo+G.ordo;
X->cote=lrg*Y.cote+G.cote;
}



void canr_loca_ferj(point *A,point *s,point *t,int N,manif_tl *loc)
{int i;
double xmi,xma,ymi,yma,zmi,zma,*h,lrg;
point *A_aux,*s_aux,*t_aux,G,temp;
homs_boun_gosm(A,3,&xmi,&xma,&ymi,&yma,&zmi,&zma);
G.absi=0.5*(xmi+xma);
G.ordo=0.5*(ymi+yma);
G.cote=0.5*(zmi+zma);
h=(double *)malloc(3*sizeof(double));
h[0]=xma-xmi;
h[1]=yma-ymi;
h[2]=zma-zmi;
lrg=0.0;
for(i=0;i<3;i++)
if(h[i]>lrg)
	lrg=h[i];
free(h);

A_aux=(point *)malloc(3*sizeof(point));
s_aux=(point *)malloc(N*sizeof(point));
t_aux=(point *)malloc(N*sizeof(point));
for(i=0;i<3;i++)
	vusz_tran_gald(A[i],G,lrg,&A_aux[i]);
for(i=0;i<N;i++)
	{vusz_tran_gald(s[i],G,lrg,&s_aux[i]);
	vusz_tran_gald(t[i],G,lrg,&t_aux[i]);
	}

wonm_loca_kind(A_aux,s_aux,t_aux,N,loc);
for(i=0;i<loc->n_grs;i++)
	{jewc_back_wosp(loc->knot[i],G,lrg,&temp);
	getf_find_rogc_todj(temp,&loc->knot[i]);
	}

free(A_aux);
free(t_aux);
free(s_aux);
}


int pizm_segm_natv(point A,point B,point X)
{int res;
double a,b,c,eps=0.0001,diff;
a=wodt_dist_gilq(A,X);
b=wodt_dist_gilq(B,X);
c=wodt_dist_gilq(A,B);
diff=fabs(a+b-c);
res=0;
if(diff<eps)
	res=1;
return res;
}


double motr_conv_rukv(point *A,point X)
{int i,nx,ts,qr,ind;
double lambda,L[3],LG,lg;
LG=0.0;   lg=0.0;
qr=1;	ind=1;
for(i=0;i<3;i++)
	{nx=i+1;
	if(nx==3)
		nx=0;
	L[i]=wodt_dist_gilq(A[i],A[nx]);
	ts=pizm_segm_natv(A[i],A[nx],X);
	if(ts==1)
		{lg=lg+wodt_dist_gilq(A[i],X);
		ind=2;
		qr=2;
		}
	if(qr==1)	
		lg=lg+L[i];
	LG=LG+L[i];
	}


if(ind==1)
	{fprintf(tmpout,"X does not belong to any edge\n");
	
	}
lambda=lg/LG;
return lambda;
}


int nebw_find_tafk(point *A,point *s,point *t,int N)
{int q=-1,i,j,r;
double lg,sm,lambda1,lambda2,lm1,lm2;
point *TEMP;

TEMP=(point *)malloc(3*sizeof(point));
for(i=0;i<3;i++)  
	{for(j=0;j<3;j++)
		{r=i+j;	if(r>=3)	r=r-3;
		getf_find_rogc_todj(A[r],&TEMP[j]);
		}
	lg=-0.1;	sm=1.1;
	for(j=0;j<N;j++)
		{lambda1=motr_conv_rukv(TEMP,s[j]);
		lambda2 =motr_conv_rukv(TEMP,t[j]);
		if((lambda1<sm)&&(lambda2>lg))
			{sm=lambda1;
			lg=lambda2;
			}
		
		if((lambda2<sm)&&(lambda1>lg))
			{sm=lambda2;
			lg=lambda1;
			}
		}
	lm1=motr_conv_rukv(TEMP,TEMP[1]);
	lm2=motr_conv_rukv(TEMP,TEMP[2]);
	if((0.0<sm)&&(sm<=lm1))
	if((0.0<lg)&&(lg>=lm2))
		{q=i;
		
		break;
		}
	}
free(TEMP);
return q;
}


