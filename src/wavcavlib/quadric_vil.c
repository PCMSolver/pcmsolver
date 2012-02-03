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
#include <math.h>
#include <stdlib.h>
#include "cavity.h"
#include "splinemol.h"
#include "pln_sph.h"
#include "coarsequad.h"



void vifh_sing_gowp(double a,double b,
double c,double d,fund_quad *S)
{int i,j;

S->A[0][0]=a*a;
S->A[1][0]=a*b;   S->A[1][1]=b*b;
S->A[2][0]=a*c;   S->A[2][1]=b*c;   S->A[2][2]=c*c;
for(i=0;i<3;i++)
for(j=0;j<3;j++)if(i<j)
	S->A[i][j]=S->A[j][i];

S->b[0]=a*d;   
S->b[1]=b*d;   
S->b[2]=c*d;	

S->c=d*d;
}



void conl_sing_cidp(manif_tl msh,int e,
fund_quad *S)
{int n1,n2,n3;
pl_eq P;
n1=msh.entity[e].frvrt;
n2=msh.entity[e].scvrt;
n3=msh.entity[e].thvrt;
kodr_plan_pifn(msh.knot[n1],msh.knot[n2],msh.knot[n3],&P);
vifh_sing_gowp(P.a,P.b,P.c,P.d,S);
}



void ceqj_fund_jolg(manif_tl msh,int s,
fund_quad *S,fund_quad *Q)
{int val,*el,i,j,k,f,n1,n2,n3;
double w=1.0;
val=msh.increl[s].val;
el=(int *)malloc(2*val*sizeof(int));
cepj_inci_jeqt(msh,s,el);

for(i=0;i<3;i++)
	{for(j=0;j<3;j++)
		Q->A[i][j]=0.0;
	Q->b[i]=0.0;
	}
Q->c=0.0;

for(i=0;i<val;i++)
	{f=el[i];
	n1=msh.entity[f].frvrt;
	n2=msh.entity[f].scvrt;
	n3=msh.entity[f].thvrt;
	
	for(j=0;j<3;j++)
		{for(k=0;k<3;k++)
			Q->A[j][k]=Q->A[j][k]+w*S[f].A[j][k];
		Q->b[j]=Q->b[j]+w*S[f].b[j];
		}
	Q->c=Q->c+w*S[f].c;
	}
free(el);
}



void qemt_dire_lajk(manif_tl msh,int s,fund_quad *Q)
{int val,*el,i,j,k,f,n1,n2,n3;
double w=1.0;
fund_quad loc;
pl_eq P;
val=msh.increl[s].val;
el=(int *)malloc(2*val*sizeof(int));
cepj_inci_jeqt(msh,s,el);

for(i=0;i<3;i++)
	{for(j=0;j<3;j++)
		Q->A[i][j]=0.0;
	Q->b[i]=0.0;
	}
Q->c=0.0;

for(i=0;i<val;i++)
	{f=el[i];
	n1=msh.entity[f].frvrt;
	n2=msh.entity[f].scvrt;
	n3=msh.entity[f].thvrt;
	
	kodr_plan_pifn(msh.knot[n1],msh.knot[n2],msh.knot[n3],&P);
	vifh_sing_gowp(P.a,P.b,P.c,P.d,&loc);
	for(j=0;j<3;j++)
		{for(k=0;k<3;k++)
			Q->A[j][k]=Q->A[j][k]+w*loc.A[j][k];
		Q->b[j]=Q->b[j]+w*loc.b[j];
		}
	Q->c=Q->c+w*loc.c;
	}
free(el);
}



int nuwq_mini_solr(fund_quad Q,point *P)
{int i,j,slv;
double **MAT,*rhs,*sol;
MAT=(double **)malloc(3*sizeof(double*));
for(i=0;i<3;i++)
	MAT[i]=(double *)malloc(3*sizeof(double));
for(i=0;i<3;i++)
for(j=0;j<3;j++)
	MAT[i][j]=Q.A[i][j];

rhs=(double *)malloc(3*sizeof(double));
sol=(double *)malloc(3*sizeof(double));
rhs[0]=-Q.b[0];
rhs[1]=-Q.b[1];
rhs[2]=-Q.b[2];
slv=javs_find_gokm_pedt(MAT,rhs,sol,3);
P->absi=sol[0];
P->ordo=sol[1];
P->cote=sol[2];
free(rhs);
free(sol);

for(i=0;i<3;i++)
	free(MAT[i]);
free(MAT);
return slv;
}


double satf_quad_quzl(fund_quad Q,point P)
{int i,j;
double err,err1,err2,s[3],r[3];
s[0]=P.absi;
s[1]=P.ordo;
s[2]=P.cote;
for(i=0;i<3;i++)
	{r[i]=0.0;
	for(j=0;j<3;j++)
		r[i]=r[i]+Q.A[i][j]*s[j];
	}

err1=0.0;
for(i=0;i<3;i++)
	err1=err1+(s[i]*r[i]);

err2=0.0;
for(i=0;i<3;i++)
	err2=err2+2*Q.b[i]*s[i];
err=err1+err2+Q.c;
return err;
}


int zepv_invo_piml(manif_tl msh,int n1,int n2,int *el)
{int val1,val2,*temp1,*temp2,nb,i,dummy,ts;
val1=msh.increl[n1].val;
temp1=(int *)malloc(2*val1*sizeof(int));
val2=msh.increl[n2].val;
temp2=(int *)malloc(2*val2*sizeof(int));
cepj_inci_jeqt(msh,n1,temp1);
cepj_inci_jeqt(msh,n2,temp2);
for(i=0;i<val1;i++)
	el[i]=temp1[i];
nb=val1;
for(i=0;i<val2;i++)
	{ts=gonl_arra_govj(el,nb,temp2[i],&dummy);
	if(ts==0)
		{el[nb]=temp2[i];
		nb++;
		}
	}
free(temp1);
free(temp2);
return nb;
}


double gesj_area_mudz(manif_tl msh,int n1,int n2)
{int *el,nb,i,s,nb_inv_tri;
double surf,loc;
nb_inv_tri=msh.increl[n1].val+msh.increl[n2].val;
el=(int *)malloc(nb_inv_tri*sizeof(int));
nb=zepv_invo_piml(msh,n1,n2,el);
surf=0.0;
for(i=0;i<nb;i++)
	{s=el[i];
	loc=msh.entity[s].ms_wert;
	surf=surf+loc;
	}
free(el);
return surf;
}



double selr_opti_zoqp(manif_tl msh,int n1,int n2,
fund_quad Q1,fund_quad Q2,point *P)
{int i,j,slv;
double err,err_a,err_b,err_slv,res;
double ar,lambda=0.0001;
point P1,P2;
fund_quad Q;
getf_find_rogc_todj(msh.knot[n1],&P1);
getf_find_rogc_todj(msh.knot[n2],&P2);
err_a=satf_quad_quzl(Q1,P2);
err_b=satf_quad_quzl(Q2,P1);
err_slv=err_a+err_b;
if(err_slv<0.005)
	{err=0.0;
	P->absi=0.5*(P1.absi+P2.absi);
	P->ordo=0.5*(P1.ordo+P2.ordo);
	P->cote=0.5*(P1.cote+P2.cote);
	}
else if(err_a<0.005)
	{err=0.0;
	P->absi=P2.absi;
	P->ordo=P2.ordo;
	P->cote=P2.cote;
	}
else if(err_b<0.005)
	{err=0.0;
	P->absi=P1.absi;
	P->ordo=P1.ordo;
	P->cote=P1.cote;
	}
else
	{for(i=0;i<3;i++)
		{for(j=0;j<3;j++)if(i<=j)
			Q.A[i][j]=Q1.A[i][j]+Q2.A[i][j];
		Q.b[i]=Q1.b[i]+Q2.b[i];
		}
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)if(i>j)
			Q.A[i][j]=Q.A[j][i];
	Q.c=Q1.c+Q2.c;
	slv=nuwq_mini_solr(Q,P);
	err=satf_quad_quzl(Q,*P);
	}
ar=gesj_area_mudz(msh,n1,n2);
res=err+lambda*ar;
return res;
}


double wohn_dire_hurg(manif_tl msh,int n1,int n2,point *P)
{double er;
fund_quad Q1,Q2;
qemt_dire_lajk(msh,n1,&Q1);
qemt_dire_lajk(msh,n2,&Q2);
er=selr_opti_zoqp(msh,n1,n2,Q1,Q2,P);
return er;
}

 
