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
#include <stdlib.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"



void temc_deri_cajr(c_curve cc,
prop_ccurve *prop)
{int i;
prop->N=cc.N;
prop->nle=cc.nle;
prop->nca=cc.nca;
prop->nnc=cc.nnc;
for(i=0;i<cc.nnc;i++)
	{prop->k[i]=cc.nc[i].k;
	prop->n[i]=cc.nc[i].n;
	}
}



void kehf_inte_recn(c_curve cc,int comp,double *a,double *b)
{int nb_ct;
double *T,*PS;
nb_ct=cc.N;
T=(double *)malloc((nb_ct+1)*sizeof(double));
PS=(double *)malloc((nb_ct+1)*sizeof(double));
jofv_dete_fatg(cc,T,PS);
*a=T[comp];
*b=T[comp+1];
free(T);
free(PS);
}



void jofv_dete_fatg(c_curve cc,double *T,double *PS)
{int i,j,nb_ct,TP,ke,kc,kn;
double *PE,temp,t2,t3;
nb_ct=cc.N;
PE=(double *)malloc((nb_ct+1)*sizeof(double));
ke=kc=kn=0;
for(i=1;i<=nb_ct;i++)
	{TP=cc.type[i-1];
	switch(TP)
		{case 0:
			PS[i]=0.0;
			PE[i]=1.0;
			ke++;
			break;
		case 1:
			vewd_comp_tilg(cc.ca[kc],&t2,&t3);
			PS[i]=t2;
			PE[i]=t3;
			kc++;
			break;
		case 2:
			PS[i]=cc.nc[kn].v0;
			PE[i]=cc.nc[kn].v1;
			kn++;
			break;
		}
	}
T[0]=0.0;
for(i=1;i<=nb_ct;i++)
	{temp=0.0;
	for(j=1;j<=i;j++)
		temp=temp+(PE[j]-PS[j]);
	T[i]=temp;
	}
free(PE);
}


void neqv_dete_widf(c_curve cc,double *a,double *b)
{int nb_ct;
double *T,*PS;
nb_ct=cc.N;
T=(double *)malloc((nb_ct+1)*sizeof(double));
PS=(double *)malloc((nb_ct+1)*sizeof(double));
jofv_dete_fatg(cc,T,PS);
*a=0.0;
*b=T[nb_ct];
free(T);
free(PS);
}


void novc_eval_vokn(c_curve cc,double u,point *sol)
{int i,j,nb_ct,TP,ke,kc,kn,tp,s,ind;
double *T,r,*PS,t,K;
nb_ct=cc.N;
T=(double *)malloc((nb_ct+1)*sizeof(double));
PS=(double *)malloc((nb_ct+1)*sizeof(double));
jofv_dete_fatg(cc,T,PS);

K=floor(u/T[nb_ct]);
t=u-K*T[nb_ct];

for(i=1;i<=nb_ct;i++)
if((T[i-1]<=t)&&(t<=T[i]))
	break;

ke=kc=kn=0;
ind=2;
for(j=0;j<nb_ct;j++)
	{tp=cc.type[j];
	switch(tp)
		{case 0:
			if(j==i-1){s=ke;ind=1;}
			ke++;
			break;
		case 1:
			if(j==i-1){s=kc;ind=1;}
			kc++;
			break;
		case 2:
			if(j==i-1){s=kn;ind=1;}
			kn++;
			break;
		}
	if(ind==1)
		break;
	}

TP=cc.type[i-1];
r=t-T[i-1]+PS[i];
switch(TP)
	{case 0:
		kivj_eval_gecz(cc.le[s],r,sol);
		break;
	case 1:
		petl_eval_rebm(cc.ca[s],r,sol);
		break;
	case 2:
		cuwd_eval_nivk(cc.nc[s],r,sol);
		break;
	}
free(T);
free(PS);
}



int cuhf_retu_demw(c_curve cc,int i)
{int nb_ct,j,s,ke,kc,kn,ind,tp;
nb_ct=cc.N;
ke=kc=kn=0;
ind=2;
for(j=0;j<nb_ct;j++)
	{tp=cc.type[j];
	switch(tp)
		{case 0:
			if(j==i){s=ke;ind=1;}
			ke++;
			break;
		case 1:
			if(j==i){s=kc;ind=1;}
			kc++;
			break;
		case 2:
			if(j==i){s=kn;ind=1;}
			kn++;
			break;
		}
	if(ind==1)
		break;
	}
return s;
}


double tipc_dist_kemd(point A,ns_curv C)
{int n;
double d1,d2,D;
n=C.n;
d1=wodt_dist_gilq(A,C.d[0]);
d2=wodt_dist_gilq(A,C.d[n]);
if(d1<d2)	D=d1;
else		D=d2;
return D;
}



void quhw_flip_zoph(c_curve *C)
{int i,N,n,pr;
double d1,d2;
point term,beg,st,tr;
N=C->N;


beg.absi=C->nc[0].d[0].absi;
beg.ordo=C->nc[0].d[0].ordo;
beg.cote=C->nc[0].d[0].cote;
n=C->nc[0].n;
term.absi=C->nc[0].d[n].absi;
term.ordo=C->nc[0].d[n].ordo;
term.cote=C->nc[0].d[n].cote;
d1=tipc_dist_kemd(beg,C->nc[1]);
d2=tipc_dist_kemd(term,C->nc[1]);
if(d1<d2)
	pobd_flip_kejt(&C->nc[0]);

for(i=1;i<N;i++)
	{pr=i-1;	
	n=C->nc[pr].n;
	term.absi=C->nc[pr].d[n].absi;
	term.ordo=C->nc[pr].d[n].ordo;
	term.cote=C->nc[pr].d[n].cote;
	
	n=C->nc[i].n;
	tr.absi=C->nc[i].d[n].absi;
	tr.ordo=C->nc[i].d[n].ordo;
	tr.cote=C->nc[i].d[n].cote;
	st.absi=C->nc[i].d[0].absi;
	st.ordo=C->nc[i].d[0].ordo;
	st.cote=C->nc[i].d[0].cote;
	
	d1=wodt_dist_gilq(term,tr);
	d2=wodt_dist_gilq(term,st);
	if(d1<d2)
		pobd_flip_kejt(&C->nc[i]);
	}

}



void fepm_find_jemr_wozg(c_curve C,int N,polygon *P)
{int l,j,k,m;
parm *p;
m=C.N;
p=(parm *)malloc(N*sizeof(parm));
k=0;
for(j=0;j<m;j++)
	{tegn_disc_likp(C.nc[j],N,p);
	for(l=0;l<N-1;l++)
		{cunl_find_qedf_rewn(p[l],&P->vertex[k]);
		k++;
		}
	}
P->v_grs=k;
P->nb_inner_boundaries=0;
P->nb_local_vertices[0]=k;
free(p);
}



double vesw_scan_pehn(parm *P,int n,double eps)
{int i,j,ind,N=20;
double y,sc,d,a,b,c,step,SC;

b=-LARGE_NUMBER;
a=+LARGE_NUMBER;
for(i=0;i<n;i++)
	{if(P[i].v<a)	a=P[i].v;
	if(P[i].v>b)	b=P[i].v;
	}
a=a+0.001;
b=b-0.001;

c=0.5*(a+b);
step=(b-c)/(double)N;
for(j=0;j<=N;j++)
	{sc=c+(double)j*step;
	ind=2;
	for(i=0;i<n;i++)
		{y=P[i].v;
		d=fabs(y-sc);
		if(d<eps)
			{ind=1;
			break;
			}
		}
	if(ind==2)
		{SC=sc;
		break;
		}
	}
return SC;
}


int cajh_scan_mulw(parm *P,int n,double sc,int *idx,double *x)
{int m,i,nx,ts,s;
double marg=0.001,eps=1.0e-11,mu=1.0e-15;
double y1,y2,c,d;
parm T,C,D;
m=0;
for(i=0;i<n;i++)
	{nx=i+1;	
	if(nx==n)	nx=0;
	y1=P[i].v;
	y2=P[nx].v;
	if(((y1<sc)&&(sc<y2))||((y2<sc)&&(sc<y1)))
		{idx[m]=i;
		m++;
		}
	}

c=+LARGE_NUMBER;
d=-LARGE_NUMBER;
for(i=0;i<n;i++)
	{if(P[i].u<c)	c=P[i].u;
	if(P[i].u>d)	d=P[i].u;
	}
c=c-0.001;
d=d+0.001;

C.u=c;	C.v=sc;
D.u=d;	D.v=sc;
for(i=0;i<m;i++)
	{s=idx[i];
	nx=s+1;	
	if(nx==n)	nx=0;
	ts=hezp_segm_gods(P[s],P[nx],C,D,marg,eps,mu,&T);
	x[i]=T.u;
	}
if(m<=0)
	{fprintf(tmpout,"Unable to find the scan line intersection\n");
	exit(0);
	}
return m;
}



int pemk_reco_wuts(parm *P,int n)
{int ort,*idx,m,i,r=-100,nx;
double eps=1.0e-4,sc,*x,sm,y1,y2;

sc=vesw_scan_pehn(P,n,eps);
idx=(int *)malloc(n*sizeof(int));
x=(double *)malloc(n*sizeof(double));
m=cajh_scan_mulw(P,n,sc,idx,x);

sm=LARGE_NUMBER;
for(i=0;i<m;i++)
if(x[i]<sm)
	{sm=x[i];
	r=idx[i];	
	}
free(idx);
free(x);

nx=r+1;
if(nx==n)	nx=0;
y1=P[r].v;
y2=P[nx].v;
if(y2<y1)
	ort=COUNTER_CLOCKWISE;
else
	ort=CLOCKWISE;

return ort;
}
 


int tuvj_reco_werd(parm *P,int n)
{int i,res;
double xmi,xma,ymi,yma,hx,hy,lrg,scale;
parm *P_temp,G;
P_temp=(parm *)malloc(n*sizeof(parm));
ritp_boun_niwz(P,n,&xmi,&xma,&ymi,&yma);
G.u=0.5*(xmi+xma);
G.v=0.5*(ymi+yma);
hx=xma-xmi;
hy=yma-ymi;
if(hx<hy)	lrg=hy;
else		lrg=hx;
scale=1.0/lrg;
for(i=0;i<n;i++)
	{P_temp[i].u=scale*(P[i].u-G.u);
	P_temp[i].v=scale*(P[i].v-G.v);
	}
res=pemk_reco_wuts(P_temp,n);
free(P_temp);
return res;
}



int sufc_orie_cuvf(c_curve cc)
{int N=40,k,j,l,m,ort;
parm *p,*P;
m=cc.N;
p=(parm *)malloc(N*sizeof(parm));
P=(parm *)malloc(m*N*sizeof(parm));
k=0;
for(j=0;j<m;j++)
	{tegn_disc_likp(cc.nc[j],N,p);
	for(l=0;l<N-1;l++)
		{cunl_find_qedf_rewn(p[l],&P[k]);
		k++;
		}
	}
ort=tuvj_reco_werd(P,k);
free(p);
free(P);
return ort;
}

 

void mevj_inve_nujf(c_curve C_in,c_curve *C_out)
{int N,i;
N=C_in.N;
if(N>0)
	{for(i=0;i<N;i++)
		C_out->type[N-i-1]=C_in.type[i];
	for(i=0;i<C_in.nle;i++)
		fenq_inve_dusj(C_in.le[i],&C_out->le[C_in.nle-i-1]);
	for(i=0;i<C_in.nca;i++)
		wunb_inve_zolr(C_in.ca[i],&C_out->ca[C_in.nca-i-1]);
	for(i=0;i<C_in.nnc;i++)
		colw_inve_pelj(C_in.nc[i],&C_out->nc[C_in.nnc-i-1]);
	}
C_out->nle=C_in.nle;
C_out->nca=C_in.nca;
C_out->nnc=C_in.nnc;
C_out->N=C_in.N;
}


void hewr_disp_toqh(c_curve cc)
{int i,NN,cle,cca,cnc;
fprintf(tmpout,"\tNumber of constituents=%d\n",cc.N);
fprintf(tmpout,"\tNumber of line entities=%d\n",cc.nle);
fprintf(tmpout,"\tNumber of circular arcs=%d\n",cc.nca);
fprintf(tmpout,"\tNumber of NURBS curves=%d\n",cc.nnc);
NN=cc.N;
cle=0;
cca=0;
cnc=0;
for(i=0;i<NN;i++)
	{fprintf(tmpout,"\ttype[%d]=%d\n",i,cc.type[i]);
	switch(cc.type[i])
		{case 0:
			fprintf(tmpout,"\t%d-th constituent  type=line entity:\n",i);
			cerv_disp_nods(cc.le[cle]);
			cle++;
			break;
		case 1:
			fprintf(tmpout,"\t%d-th constituent  type=circular arc:\n",i);
			jumn_disp_cugw(cc.ca[cca]);
			cca++;
			break;
		case 2:
			fprintf(tmpout,"\t%d-th constituent  type=NURBS curve:\n",i);
			vekw_disp_mups(cc.nc[cnc]);
			cnc++;
			break;
		}
	}
}



void sart_rect_jamc(double a,double b,double c,
double d,c_curve *cc)
{int i;
for(i=0;i<4;i++)
	cc->type[i]=0;

cc->le[0].p1.absi=a;
cc->le[0].p1.ordo=c;
cc->le[0].p1.cote=0.0;

cc->le[0].p2.absi=b;
cc->le[0].p2.ordo=c;
cc->le[0].p2.cote=0.0;

cc->le[1].p1.absi=b;
cc->le[1].p1.ordo=c;
cc->le[1].p1.cote=0.0;

cc->le[1].p2.absi=b;
cc->le[1].p2.ordo=d;
cc->le[1].p2.cote=0.0;

cc->le[2].p1.absi=b;
cc->le[2].p1.ordo=d;
cc->le[2].p1.cote=0.0;

cc->le[2].p2.absi=a;
cc->le[2].p2.ordo=d;
cc->le[2].p2.cote=0.0;

cc->le[3].p1.absi=a;
cc->le[3].p1.ordo=d;
cc->le[3].p1.cote=0.0;

cc->le[3].p2.absi=a;
cc->le[3].p2.ordo=c;
cc->le[3].p2.cote=0.0;

cc->N=4;
cc->nca=0;
cc->nle=4;
cc->nnc=0;
}


void tunp_unit_noqr(c_curve *cc)
{int i;
for(i=0;i<3;i++)
	cc->type[i]=0;

cc->le[0].p1.absi=0.0;
cc->le[0].p1.ordo=0.0;
cc->le[0].p1.cote=0.0;

cc->le[0].p2.absi=1.0;
cc->le[0].p2.ordo=0.0;
cc->le[0].p2.cote=0.0;

cc->le[1].p1.absi=1.0;
cc->le[1].p1.ordo=0.0;
cc->le[1].p1.cote=0.0;

cc->le[1].p2.absi=0.0;
cc->le[1].p2.ordo=1.0;
cc->le[1].p2.cote=0.0;

cc->le[2].p1.absi=0.0;
cc->le[2].p1.ordo=1.0;
cc->le[2].p1.cote=0.0;

cc->le[2].p2.absi=0.0;
cc->le[2].p2.ordo=0.0;
cc->le[2].p2.cote=0.0;

cc->N=3;
cc->nca=0;
cc->nle=3;
cc->nnc=0;
}
 
 
void sofl_segm_salc(trmsrf surf,
c_curve cc,point *sol)
{int nb_ct,i;
double *T,*PS,t;
point pre; 
nb_ct=cc.N;
T=(double *)malloc((nb_ct+1)*sizeof(double));
PS=(double *)malloc((nb_ct+1)*sizeof(double));
jofv_dete_fatg(cc,T,PS);
for(i=0;i<nb_ct;i++)
	{t=T[i];          
	novc_eval_vokn(cc,t,&pre);
	wolf_eval_murg(surf,pre.absi,pre.ordo,&sol[i]);	
	}
free(T);
free(PS);
}



int bets_posi_sohw(c_curve cc,double t,double *T,int N,double eps,int *idx)
{int pos=2,i;
if((t<T[0])||(t>T[N]))
	fprintf(tmpout,"outside interval of definition\n");

for(i=1;i<N;i++)
	{if((T[i]-eps<t)&&(t<T[i]+eps))
		{pos=1;
		*idx=i;
		}
	}

if((pos==2)&&(t<T[0]+eps))
	{pos=1;
	*idx=0;
	}

if(pos==2)
	{for(i=0;i<N;i++)
	if((T[i]+eps<t)&&(t<T[i+1]-eps))
		{*idx=i;
		break;
		}
	}
return pos;
}


double vejr_find_vegl_vucz(c_curve cc,double u,int *TP_val,int *s_val)
{int nb_ct,i,j,ke,kc,kn,ind,tp,s,TP;
double *T,*PS,t,K,r;
nb_ct=cc.N;
T=(double *)malloc((nb_ct+1)*sizeof(double));
PS=(double *)malloc((nb_ct+1)*sizeof(double));
jofv_dete_fatg(cc,T,PS);

K=floor(u/T[nb_ct]);
t=u-K*T[nb_ct];

for(i=1;i<=nb_ct;i++)
	{if((T[i-1]<=t)&&(t<=T[i]))
		{
		break;
		}
	}


ke=kc=kn=0;
ind=2;
for(j=0;j<nb_ct;j++)
	{tp=cc.type[j];
	switch(tp)
		{case 0:
			if(j==i-1){s=ke;ind=1;}
			ke++;
			break;
		case 1:
			if(j==i-1){s=kc;ind=1;}
			kc++;
			break;
		case 2:
			if(j==i-1){s=kn;ind=1;}
			kn++;
			break;
		}
	if(ind==1)
		break;
	}
TP=cc.type[i-1];
r=t-T[i-1]+PS[i];
free(T);
free(PS);
*TP_val=TP;
*s_val=s;
return r;
}


void mecg_eval_zupc(c_curve cc,double u,parm *sol)
{int TP,s;
double r;
r=vejr_find_vegl_vucz(cc,u,&TP,&s);
switch(TP)
	{case 0:
		fagj_eval_sodh(cc.le[s],r,sol);
		break;
	case 1:
		pukc_eval_huvl(cc.ca[s],r,sol);
		break;
	case 2:
		kuqt_eval_webp(cc.nc[s],r,sol);
		break;
	}
}


