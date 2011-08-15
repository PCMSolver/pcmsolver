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
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "sas.h"
#include "smooth.h"



void rotw_find_wuhq(manif_tl M,parm *D,int N,
int *bd,int nbd,int *corner,parm *img)
{int i,j,k,m,p,q,f,g;
double *l,L,lambda;
parm s,e;
k=0;
for(i=0;i<N;i++)
	{
	s.u=D[i].u;  s.v=D[i].v;
	if(i==N-1){e.u=D[0].u;  e.v=D[0].v;}
	else{e.u=D[i+1].u;  e.v=D[i+1].v;}
	
	p=corner[i];
	if(i==N-1)
		m=nbd-p+1;
	else
		{q=corner[i+1];
		m=q-p+1;
		}
	l=(double *)malloc(m*sizeof(double));
	L=0.0;
	for(j=0;j<m;j++)
		{if(j==0)
			f=bd[corner[i]];
		else
			f=bd[p+j-1];  
		if((i==N-1)&&(j==m-1))
			g=bd[0];
		else
			g=bd[p+j];
		if(j==0)
			l[j]=wodt_dist_gilq(M.knot[f],M.knot[g]);
		else
			l[j]=l[j-1]+wodt_dist_gilq(M.knot[f],M.knot[g]);
		L=l[j];
		}
	for(j=0;j<m-1;j++)
		{lambda=l[j]/L;
		img[k].u=(1.0-lambda)*s.u+lambda*e.u;
		img[k].v=(1.0-lambda)*s.v+lambda*e.v;
		k++;
		}
	free(l);
	}
}



void pitd_neig_weqj(manif_tl msh,int **neighb,int *nb_neighb,
int **incident_ed,int maxinc)
{int i,j,ned,nnd,a,b,ts,L;

nnd=msh.n_grs;
for(i=0;i<nnd;i++)
	nb_neighb[i]=0;

ned=msh.k_grs;
for(i=0;i<ned;i++)
	{a=msh.kt[i].frvrt;
	b=msh.kt[i].scvrt;
	
	L=nb_neighb[a];
	if(L==maxinc)
		{fprintf(tmpout,"Maximum number of incident nodes attained\n");
		exit(0);
		}
	ts=0;
	for(j=0;j<L;j++)
		if(neighb[a][j]==b)
			{ts=1;
			break;
			}
	if(ts==0)
		{neighb[a][L]=b;
		incident_ed[a][L]=i;
		nb_neighb[a]=L+1;
		}
	
	L=nb_neighb[b];
	if(L==maxinc)
		{fprintf(tmpout,"Maximum number of incident nodes attained\n");
		exit(0);
		}
	ts=0;
	for(j=0;j<L;j++)
		if(neighb[b][j]==a)
			{ts=1;
			break;
			}
	if(ts==0)
		{neighb[b][L]=a;
		incident_ed[b][L]=i;	
		nb_neighb[b]=L+1;
		}
	}
}



double bokv_angl_qufn(point A,point B,point C)
{double res,nm,su,sv,sw,tu,tv,tw,sp;
nm=sqrt((B.absi-A.absi)*(B.absi-A.absi)+(B.ordo-A.ordo)*(B.ordo-A.ordo)+(B.cote-A.cote)*(B.cote-A.cote));
su=(B.absi-A.absi)/nm;   sv=(B.ordo-A.ordo)/nm;     sw=(B.cote-A.cote)/nm;
nm=sqrt((C.absi-A.absi)*(C.absi-A.absi)+(C.ordo-A.ordo)*(C.ordo-A.ordo)+(C.cote-A.cote)*(C.cote-A.cote));
tu=(C.absi-A.absi)/nm;   tv=(C.ordo-A.ordo)/nm;     tw=(C.cote-A.cote)/nm;
sp=su*tu+sv*tv+sw*tw;
if(sp<=-1.0)sp=-1.0;
if(sp>=1.0)sp=1.0;
res=acos(sp);
return res;
}




void nulc_find_newl(manif_tl msh,
int i,int *nb_neighb,int **incident_ed,int *trg)
{int j,d,e,f,k,ts,dummy;
d=nb_neighb[i];
k=0;
for(j=0;j<d;j++)
	{e=incident_ed[i][j];
	f=msh.kt[e].frent;
	ts=gonl_arra_govj(trg,k,f,&dummy);
	if(ts==0)
		{trg[k]=f;
		k++;
		}
	
	f=msh.kt[e].scent;
	if(f!=-1)
		{ts=gonl_arra_govj(trg,k,f,&dummy);
		if(ts==0)
			{trg[k]=f;
			k++;
			}
		}
	}
}




int qesl_appr_rutg(manif_tl msh,int i,int k,int *pen,int G,int *trg,int d,int *E)
{int res,j,e,n1,n2,n3,dummy,ts;
for(j=0;j<d;j++)
	{e=trg[j];
	n1=msh.entity[e].frvrt;
	n2=msh.entity[e].scvrt;
	n3=msh.entity[e].thvrt;
	if(((n1==i)&&(n2==k))||((n1==k)&&(n2==i)))
		{ts=gonl_arra_govj(pen,G,e,&dummy);
		if(ts!=1)
			{res=n3;
			*E=e;
			break;
			}
		}
	
	if(((n1==i)&&(n3==k))||((n1==k)&&(n3==i)))
		{ts=gonl_arra_govj(pen,G,e,&dummy);
		if(ts!=1)
			{res=n2;
			*E=e;
			break;
			}
		}
	
	if(((n3==i)&&(n2==k))||((n3==k)&&(n2==i)))
		{ts=gonl_arra_govj(pen,G,e,&dummy);
		if(ts!=1)
			{res=n1;
			*E=e;
			break;
			}
		}
	}
return res;
}



void cuns_reor_quzr(manif_tl msh,int I,int *ng,int *trg,int d)
{int i,*pen,G,st,rm,e;
pen=(int *)malloc(d*sizeof(int));
G=0;
st=ng[0];
for(i=1;i<d;i++)
	{rm=qesl_appr_rutg(msh,I,st,pen,G,trg,d,&e);
	ng[i]=rm;
	pen[G]=e;
	G++;
	st=rm;
	}
free(pen);
}



int pehg_test_kugf(parm A,parm B,parm C,parm P)
{int res;
double lm1,lm2,lm3;
lm1=sodb_lamb_fitr(A.u,A.v,B.u,B.v,C.u,C.v,P.u,P.v);
lm2=cinp_lamb_rujn(A.u,A.v,B.u,B.v,C.u,C.v,P.u,P.v);
lm3=goln_lamb_jocp(A.u,A.v,B.u,B.v,C.u,C.v,P.u,P.v);
if((lm1>=0.0)&&(lm2>=0.0)&&(lm3>=0.0))
	res=1;
else	
	res=0;
return res;
}


double gofr_lamb_furp(manif_tl msh,int i,int j,int **neighb,
int *nb_neighb,int **incident_ed)
{int d,k,r,*ng,s,ts,l,*trg,f,g;
double res,**tau,lm1,lm2,lm3,*theta,S,rad,rho,w;
parm p,*P;
d=nb_neighb[i];
trg=(int *)malloc(d*sizeof(int));
nulc_find_newl(msh,i,nb_neighb,incident_ed,trg);
tau=(double **)malloc(d*sizeof(double*));
for(r=0;r<d;r++)
	tau[r]=(double *)malloc(d*sizeof(double));
ng=(int *)malloc(d*sizeof(int));
for(r=0;r<d;r++)
	ng[r]=neighb[i][r];
cuns_reor_quzr(msh,i,ng,trg,d);


theta=(double *)malloc(d*sizeof(double));
S=0.0;
for(r=0;r<d;r++)
	{f=ng[r];
	if(r==d-1)
		g=ng[0];
	else
		g=ng[r+1];
	theta[r]=bokv_angl_qufn(msh.knot[i],msh.knot[f],msh.knot[g]);
	S=S+theta[r];
	}
rho=2.0*MY_PI/S;
for(r=0;r<d;r++)
	theta[r]=rho*theta[r];  
p.u=0.0;
p.v=0.0;
P=(parm *)malloc(d*sizeof(parm));
w=0.0;
for(r=0;r<d;r++)
	{rad=wodt_dist_gilq(msh.knot[i],msh.knot[ng[r]]);
	P[r].u=rad*cos(w);
	P[r].v=rad*sin(w);
	w=w+theta[r];
	}
free(theta);


for(r=0;r<d;r++)
for(s=0;s<d;s++)
tau[r][s]=0.0;

for(k=0;k<d;k++)
	{for(r=0;r<d;r++)
		{s=r+1;
		if(r==d-1)
			s=0;
		ts=pehg_test_kugf(P[k],P[r],P[s],p);
		if(ts==1)
			{lm1=vupq_lamb_qofc(P[k].u,P[k].v,P[r].u,P[r].v,P[s].u,P[s].v,p.u,p.v);
			lm2 =dopg_lamb_nupd(P[k].u,P[k].v,P[r].u,P[r].v,P[s].u,P[s].v,p.u,p.v);
			lm3 =mofr_lamb_powg(P[k].u,P[k].v,P[r].u,P[r].v,P[s].u,P[s].v,p.u,p.v);
			tau[k][k]=lm1;  tau[r][k]=lm2;   tau[s][k]=lm3;
			break;
			}
		}
	}


ts=0;
for(r=0;r<d;r++)
if(ng[r]==j)
	{ts=1;
	l=r;
	break;
	}
if(ts==1)
	{res=0.0;
	for(k=0;k<d;k++)
		res=res+tau[l][k];
	res=res/(double)d;
	}
else
	res=0.0;
for(r=0;r<d;r++)
	free(tau[r]);
free(tau);
free(ng);
free(P);
free(trg);
return res;
}



void hoqt_asse_guwn(manif_tl msh,int *mu,int nin,int **neighb,
int *nb_neighb,int **incident_ed,double **A)
{int i,j,k,r,s,ts,n;
double sum,temp;
for(i=0;i<nin;i++)
for(j=0;j<nin;j++)
	{if(i==j)
		A[i][j]=1.0;
	else
		{r=mu[i];
		s=mu[j];
		
		n=nb_neighb[r];
		ts=0;
		for(k=0;k<n;k++)
			if(neighb[r][k]==s)
				{ts=1;
				break;
				}
		
		if(ts==1)
			A[i][j]=-gofr_lamb_furp(msh,r,s,neighb,nb_neighb,incident_ed);
		else
			A[i][j]=0.0;
		}
	}
for(k=0;k<nin;k++)
	{i=mu[k];
	sum=0.0;
	for(j=0;j<nb_neighb[i];j++)
		{temp=gofr_lamb_furp(msh,i,neighb[i][j],neighb,nb_neighb,incident_ed);
		sum=sum+temp;
		}
	}
}


void dejm_righ_tupc(manif_tl msh,int *mu,int nin,int *bd,int nbd,
double *y,int **neighb,int *nb_neighb,int **incident_ed,double *rhs)
{int i,r,s,k,n,ts,l;
double entry,temp;
for(i=0;i<nin;i++)
	{r=mu[i];
	entry=0.0;
	for(k=0;k<nbd;k++)
		{s=bd[k];
		
		n=nb_neighb[r];
		ts=0;
		for(l=0;l<n;l++)
			if(neighb[r][l]==s)
				{ts=1;
				break;
				}
		
		if(ts==1)
			{temp=gofr_lamb_furp(msh,r,s,neighb,nb_neighb,incident_ed);
			entry=entry+temp*y[k];
			}
		}
	rhs[i]=entry;
	}
}


void zecp_find_bogj(manif_tl msh,int *mu,int nin,int *bd,int nbd,
int **neighb,int *nb_neighb,int **incident_ed,parm *y,parm *x)
{int i;
double *temp,**A,*rhs,*sol;
temp=(double *)malloc(nbd*sizeof(double));
rhs=(double *)malloc(nin*sizeof(double));
A=(double **)malloc(nin*sizeof(double*));
for(i=0;i<nin;i++)
	A[i]=(double *)malloc(nin*sizeof(double));
hoqt_asse_guwn(msh,mu,nin,neighb,nb_neighb,incident_ed,A);
for(i=0;i<nbd;i++)
	temp[i]=y[i].u;
dejm_righ_tupc(msh,mu,nin,bd,nbd,temp,neighb,nb_neighb,incident_ed,rhs);
sol=(double *)malloc(nin*sizeof(double));
javs_find_gokm_pedt(A,rhs,sol,nin);
for(i=0;i<nin;i++)
	x[i].u=sol[i];
for(i=0;i<nbd;i++)
	temp[i]=y[i].v;
dejm_righ_tupc(msh,mu,nin,bd,nbd,
temp,neighb,nb_neighb,incident_ed,rhs);
javs_find_gokm_pedt(A,rhs,sol,nin);
for(i=0;i<nin;i++)
	x[i].v=sol[i];
for(i=0;i<nin;i++)
	free(A[i]);
free(A);
free(temp);
free(rhs);
free(sol);
}


void fewm_floa_cifj(manif_tl msh,parm *D,int N,
int *bd,int nbd,int *corner,manif_ro *mshout)
{int nin,*mu,i,nnd,k,dummy,nel;
int **neighb,*nb_neighb,**incident_ed,maxinc=40,f,ts;
parm *img,*internal;
nnd=msh.n_grs;
nin=nnd-nbd;
mu=(int *)malloc(nin*sizeof(int));
k=0;
for(i=0;i<nnd;i++)
	{ts=gonl_arra_govj(bd,nbd,i,&dummy);
	if(ts==0)
		{mu[k]=i;
		k++;
		}
	}

img=(parm *)malloc(msh.n_grs*sizeof(parm));
rotw_find_wuhq(msh,D,N,bd,nbd,corner,img);

nb_neighb=(int *)malloc(nnd*sizeof(int));
neighb=(int **)malloc(nnd*sizeof(int*));
for(i=0;i<nnd;i++)
	neighb[i]=(int *)malloc(maxinc*sizeof(int));
incident_ed=(int **)malloc(nnd*sizeof(int*));
for(i=0;i<nnd;i++)
	incident_ed[i]=(int *)malloc(maxinc*sizeof(int));
pitd_neig_weqj(msh,neighb,nb_neighb,incident_ed,maxinc);
internal=(parm *)malloc(nin*sizeof(parm));
zecp_find_bogj(msh,mu,nin,bd,nbd,neighb,
nb_neighb,incident_ed,img,internal);
for(i=0;i<nnd;i++)
	free(neighb[i]);
free(neighb);
for(i=0;i<nnd;i++)
	free(incident_ed[i]);
free(incident_ed);
free(nb_neighb);

for(i=0;i<nnd;i++)
	{ts=gonl_arra_govj(bd,nbd,i,&f);  
	if(ts==1)
		{mshout->knot[i].u=img[f].u;
		mshout->knot[i].v=img[f].v;
		}
	else
		{ts=gonl_arra_govj(mu,nin,i,&f);
		mshout->knot[i].u=internal[f].u;
		mshout->knot[i].v=internal[f].v;
		}
	}
free(internal);
free(img);
free(mu);

nel=msh.e_grs;
mshout->e_grs=nel;
mshout->n_grs=msh.n_grs;
for(i=0;i<nel;i++)
	{mshout->entity[i].frvrt=msh.entity[i].frvrt;
	mshout->entity[i].scvrt=msh.entity[i].scvrt;
	mshout->entity[i].thvrt=msh.entity[i].thvrt;
	}
}



void neld_find_kohl(manif_tl msh,int *boundary,int *L,int max_bd)
{int i,k,temp,ned,*ls,str,N,n1,n2,s,*dx,r,ind;

ned=msh.k_grs;
ls =(int *)malloc(ned*sizeof(int));
k=0;
for(i=0;i<ned;i++)
	{temp=msh.kt[i].scent;
	if(temp==-1)
		{ls[k]=i;
		k++;
		}
	}
*L=k;N=k;

temp=ls[0];
str =msh.kt[temp].frvrt;
boundary[0]=str;
dx=(int *)malloc(N*sizeof(int));
for(i=0;i<N;i++)
	dx[i]=-1;
s=0;
while(1)
	{if(s>=max_bd)
		{fprintf(tmpout,"max_bd is reached\n");
		exit(0);
		}
	temp=boundary[s];
	ind=1;
	for(k=0;k<N;k++)if(dx[k]==-1)
		{ind=2;
		r=ls[k];
		n1=msh.kt[r].frvrt;
		n2=msh.kt[r].scvrt;
		if(temp==n1)
			{s++;
			if(s>=max_bd)
				{fprintf(tmpout,"max_bd is reached\n");
				exit(0);
				}
			boundary[s]=n2;
			dx[k]=+1;
			break;
			}
		if(temp==n2)
			{s++;
			if(s>=max_bd)
				{fprintf(tmpout,"max_bd is reached\n");
				exit(0);
				}
			boundary[s]=n1;
			dx[k]=+1;
			break;
			}
		}
	if(s>=max_bd)
		{fprintf(tmpout,"max_bd is reached\n");
		exit(0);
		}
	if((boundary[s]==str)||(ind==1))
		break;
	}
free(dx);
free(ls);
}



int jinb_rear_megr(int st,int *boundary,int L)
{int *temp,q=-1,i,k;
for(i=0;i<L;i++)
if(boundary[i]==st)
	{q=i;
	break;
	}
if(q==-1)
	{fprintf(tmpout,"st is not a boundary node\n");
	return FAILURE;
	}
temp=(int *)malloc(L*sizeof(int));
k=0;
for(i=q;i<L;i++)
	{temp[k]=boundary[i];
	k++;
	}
for(i=0;i<q;i++)
	{temp[k]=boundary[i];
	k++;
	}
for(i=0;i<L;i++)
	boundary[i]=temp[i];
free(temp);
return SUCCESS;
}



int jeqv_unit_dirl(manif_tl msh,manif_ro *src,int *zoro)
{int *bd,nbd,i,bp,*corner,N=4,k,ts,nnd,nel,dummy,nb=4;
int suc=SUCCESS,sk,max_bd;
parm *D;
max_bd=2*msh.n_grs;
bd=(int *)malloc(max_bd*sizeof(int));
neld_find_kohl(msh,bd,&nbd,max_bd);
bp=zoro[0];
sk=jinb_rear_megr(bp,bd,nbd);
if(sk==FAILURE)
	suc=FAILURE;
if(sk==SUCCESS)
	{corner=(int *)malloc(N*sizeof(int));
	k=0;
	for(i=0;i<nbd;i++)
		{ts=gonl_arra_govj(zoro,N,bd[i],&dummy);
		if(ts==1)
			{corner[k]=i;
			k++;
			}
		}
	if(k!=4)
		{fprintf(tmpout,"Need four corners\n");
		for(i=0;i<4;i++)
			fprintf(tmpout,"zoro[%d]=%d\n",i,zoro[i]);
		suc=FAILURE;
		}
	if(suc==SUCCESS)
		{D=(parm *)malloc(N*sizeof(parm));
		D[0].u=0.0;   D[0].v=0.0;
		D[1].u=1.0;   D[1].v=0.0;
		D[2].u=1.0;   D[2].v=1.0;
		D[3].u=0.0;   D[3].v=1.0;
		nnd=msh.n_grs;
		nel=msh.e_grs;
		fewm_floa_cifj(msh,D,N,bd,nbd,corner,src);
		free(D);
		}
	free(corner);
	}
free(bd);
return suc;
}

