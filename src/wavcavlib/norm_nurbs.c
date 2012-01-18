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
#include <math.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "eval.h"
#include "smooth.h"


void fojb_allo_qugd(prop_n_curv pnc,ns_curv1D *nc)
{int nn,kk;
nn     =pnc.n;
kk     =pnc.k;
nc->w  =(double *)malloc((nn+1)*sizeof(double));
nc->d  =(double *)malloc((nn+1)*sizeof(double));
nc->tau=(double *)malloc((nn+kk+1)*sizeof(double));
nc->k=kk;
nc->n=nn;
}


void gojt_dest_wujl(ns_curv1D *nc)
{int nn,kk;
nn=nc->n;
kk=nc->k;
if(nn>0)
	{free(nc->w);
	free(nc->d);
	free(nc->tau);
	}
}


void juvm_allo_sehv(prop_n_surf pns,ns_surf *ns)
{int i,nnu,nnv,kku,kkv;
nnu=pns.nu;
nnv=pns.nv;
kku=pns.ku;
kkv=pns.kv;
ns->w=(double **)malloc((nnu+1)*sizeof(double*)); 
for(i=0;i<nnu+1;i++) 
	ns->w[i]=(double *)malloc((nnv+1)*sizeof(double)); 
ns->d=(point **)malloc((nnu+1)*sizeof(point)); 
for(i=0;i<nnu+1;i++) 
	ns->d[i]=(point *)malloc((nnv+1)*sizeof(point)); 
ns->frknt=(double *)malloc((nnu+kku+1)*sizeof(double)); 
ns->scknt=(double *)malloc((nnv+kkv+1)*sizeof(double)); 
ns->ku   =kku;
ns->kv   =kkv;
ns->nu   =nnu;
ns->nv   =nnv;
}



void destroy_nurbs_surface_alx(prop_n_surf pns,
ns_surf *ns)
{int i,nnu,nnv,kku,kkv;
nnu=pns.nu;
nnv=pns.nv;
kku=pns.ku;
kkv=pns.kv;
free(ns->frknt);
free(ns->scknt);
for(i=0;i<nnu+1;i++) 
	free(ns->d[i]);
free(ns->d);
for(i=0;i<nnu+1;i++) 
	free(ns->w[i]);
free(ns->w);
}


point ** allocate_mat_point(int n,int m)
{int i;
point **M;
M=(point **)malloc(n*sizeof(point));
for(i=0;i<n;i++)
	M[i]=(point *)malloc(m*sizeof(point));
return M;
}


void dirj_free_sukl(point **M,int n,int m)
{int i;
for(i=n-1;i>=0;i--)
	free((char *)M[i]);
free(M);
}


double jurw_eval_voch(ns_curv1D nc,double t)
{int k,n;
double res;
k=nc.k;
n=nc.n;
res=husv_curv_kogc(k,n,nc.d,nc.tau,t);




return res;
}


double vulj_surf_mukl(int ku,int kv,int nu,int nv,double **d,
double *tau_u,double *tau_v,double u,double v,int r_u,int r_v)
{int i,j;
double *temp1,*temp2,sol;
temp1=(double *)malloc((nu+1)*sizeof(double));
temp2=(double *)malloc((nv+1)*sizeof(double)); 
for(i=0;i<=nu;i++)
	{for(j=0;j<=nv;j++)
		temp2[j]=d[i][j];
	temp1[i]=zopg_curv_zecg(kv,nv,r_v,temp2,tau_v,v);
	}
sol=zopg_curv_zecg(ku,nu,r_u,temp1,tau_u,u);
free(temp1);
free(temp2);
return sol;
}


void golj_surf_lonr(int ku,int kv,int nu,int nv,point **d,
double *tau_u,double *tau_v,double u,double v,point *sol)
{int i,j,r_u,r_v;
double **temp1,**temp2,**temp3;
r_u  =pekq_knot_zevj(u,nu,ku,tau_u);
r_v  =pekq_knot_zevj(v,nv,kv,tau_v);
temp1=allocate_mat(nu+1,nv+1);
temp2=allocate_mat(nu+1,nv+1);
temp3=allocate_mat(nu+1,nv+1);
for(i=0;i<=nu;i++)
for(j=0;j<=nv;j++)
	{temp1[i][j]=d[i][j].absi;
	temp2[i][j]=d[i][j].ordo;
	temp3[i][j]=d[i][j].cote;
	}
sol->absi=vulj_surf_mukl(ku,kv,nu,nv,temp1,tau_u,tau_v,u,v,r_u,r_v);
sol->ordo=vulj_surf_mukl(ku,kv,nu,nv,temp2,tau_u,tau_v,u,v,r_u,r_v);
sol->cote=vulj_surf_mukl(ku,kv,nu,nv,temp3,tau_u,tau_v,u,v,r_u,r_v);
tehg_free_dacp(temp1,nu+1,nv+1);
tehg_free_dacp(temp2,nu+1,nv+1);
tehg_free_dacp(temp3,nu+1,nv+1);
}


double lohk_surf_gahk(int ku,int kv,int nu,int nv,double **d,
double *tau_u,double *tau_v,double u,double v)
{int r_u,r_v;
double sol;

r_u=pekq_knot_zevj(u,nu,ku,tau_u);
r_v=pekq_knot_zevj(v,nv,kv,tau_v);

sol=vulj_surf_mukl(ku,kv,nu,nv,d,tau_u,tau_v,u,v,r_u,r_v);
return sol;
}


void decr_find_fiht_pefl(int ku,int kv,int nu,int nv,point **d,double **w,
double *tau_u,double *tau_v,double u,double v,point *sol)
{int i,j;
double den;
point **D,num;
D=(point **)malloc((nu+1)*sizeof(point));
for(i=0;i<nu+1;i++)
	D[i]=(point *)malloc((nv+1)*sizeof(point));
for(i=0;i<=nu;i++)
for(j=0;j<=nv;j++)
	{D[i][j].absi=d[i][j].absi*w[i][j];
	D[i][j].ordo=d[i][j].ordo*w[i][j];
	D[i][j].cote=d[i][j].cote*w[i][j];
	}
golj_surf_lonr(ku,kv,nu,nv,D,tau_u,tau_v,u,v,&num);
den=lohk_surf_gahk(ku,kv,nu,nv,w,tau_u,tau_v,u,v);
for(i=0;i<nu+1;i++)
	free(D[i]);
free(D);
sol->absi=num.absi/den;
sol->ordo=num.ordo/den;
sol->cote=num.cote/den;
}


void cilj_eval_qelf(ns_surf ns,double u,double v,point *sol)
{int p3,kku,kkv,nnu,nnv;
p3=ns.prop3;
switch(p3)
	{case 0:
		kku=ns.ku;    nnu=ns.nu;
		kkv=ns.kv;    nnv=ns.nv;
		decr_find_fiht_pefl(kku,kkv,nnu,nnv,ns.d,ns.w,ns.frknt,ns.scknt,u,v,sol);
		break;
	case 1:
		kku=ns.ku;    nnu=ns.nu;
		kkv=ns.kv;    nnv=ns.nv;
		golj_surf_lonr(kku,kkv,nnu,nnv,ns.d,ns.frknt,ns.scknt,u,v,sol);
		break;
	}
}

 
void qucp_find_pogc_gecz(ns_surf nsi,ns_surf *nso)
{int nnu,nnv,kku,kkv,i,j;
nnu     = nsi.nu;  
nnv     = nsi.nv;  
kku     = nsi.ku;  
kkv     = nsi.kv;
nso->nu = nnu;  
nso->nv = nnv;  
nso->ku = kku;  
nso->kv = kkv;

for(i=0;i<=nnu;i++)
for(j=0;j<=nnv;j++)
	{nso->w[i][j]=nsi.w[i][j];
	getf_find_rogc_todj(nsi.d[i][j],&(nso->d[i][j]));
	}

for(i=0;i<=nnu+kku;i++)
	nso->frknt[i]=nsi.frknt[i];
for(i=0;i<=nnv+kkv;i++)
	nso->scknt[i]=nsi.scknt[i];

nso->u0=nsi.u0;  
nso->u1=nsi.u1;
nso->v0=nsi.v0;  
nso->v1=nsi.v1;

nso->prop1=nsi.prop1;
nso->prop2=nsi.prop2;
nso->prop3=nsi.prop3;
nso->prop4=nsi.prop4;
nso->prop5=nsi.prop5;
}



void carw_flip_vewk(ns_surf *S)
{int nu,nv,ku,kv,i,j;
double a,b,df;
prop_n_surf pns;
ns_surf temp;
nu=S->nu;	ku=S->ku;
nv=S->nv;	kv=S->kv;
pns.nu=nu;	pns.ku=ku;
pns.nv=nv;	pns.kv=kv;
juvm_allo_sehv(pns,&temp);
for(i=0;i<=nu;i++)
for(j=0;j<=nv;j++)
	{getf_find_rogc_todj(S->d[i][nv-j],&temp.d[i][j]);
	temp.w[i][j]=S->w[i][nv-j];
	}
a=S->v0;	b=S->v1;
for(i=0;i<=nv+kv;i++)
	{df=b-S->scknt[i];
	temp.scknt[nv+kv-i]=a+df;
	}

for(i=0;i<=nu;i++)
for(j=0;j<=nv;j++)
	{getf_find_rogc_todj(temp.d[i][j],&S->d[i][j]);
	S->w[i][j]=temp.w[i][j];
	}
for(i=0;i<=nv+kv;i++)
	S->scknt[i]=temp.scknt[i];
destroy_nurbs_surface_alx(pns,&temp);

}



int current_norm_dir(ns_surf *S,int n_patches,int m)
{int parity,p,i,j,k,n_loc,n_sam,res;
int p_spec,i_spec,k_spec,n_left,n_right;
double x,y,step,lrg,cte,sp;
point *P,img,pt_dap;
vect3D N,dir;
parm *q;

if(n_patches<=1)
	return OUTWARD_ORIENT;
parity=m % 2;
if(parity==0)
	{fprintf(tmpout,"We should have an odd number\n");
	exit(0);
	}

step=1.0/((double)m+1.0);
q=(parm *)malloc(m*m*sizeof(parm));
k=0;
for(i=0;i<m;i++)
	{x=((double)i+1.0)*step;
	for(j=0;j<m;j++)
		{y=((double)j+1.0)*step;
		q[k].u=x;
		q[k].v=y;
		k++;
		}
	}
n_loc=k;

lrg=-LARGE_NUMBER;
P=(point *)malloc(n_patches*m*m*sizeof(point));
k=0;
for(p=0;p<n_patches;p++)
for(i=0;i<n_loc;i++)
	{x=q[i].u;
	y=q[i].v;
	cilj_eval_qelf(S[p],x,y,&img);
	getf_find_rogc_todj(img,&P[k]);
	cte=P[k].absi;
	if(cte>lrg)
		{p_spec=p;
		i_spec=i;
		k_spec=k;
		lrg=cte;
		}
	k++;
	}
n_sam=k;
x=q[i_spec].u;
y=q[i_spec].v;
gurn_norm_tegm(S[p_spec],x,y,1.0e-4,&N);
getf_find_rogc_todj(P[k_spec],&pt_dap);
free(q);

n_left=0;
n_right=0;
for(i=0;i<n_sam;i++)
if(i!=k_spec)
	{bofp_form_nukv(pt_dap,P[i],&dir);
	sp=rocv_scal_toqc(dir,N);
	if(sp>0.0)	n_right++;
	else		n_left++;
	}
free(P);
fprintf(tmpout,"n_left=%d   n_right=%d\n",n_left,n_right);
if(n_left>n_right)	
	res=OUTWARD_ORIENT;
else
	res=INWARD_ORIENT;
return res;
}

  
void automatic_orientation(ns_surf *S,int n_patches,int desired_ort)
{int cur,m=3,i;
cur=current_norm_dir(S,n_patches,m);


if(cur!=desired_ort)
	{for(i=0;i<n_patches;i++)
		carw_flip_vewk(&S[i]);
	}
}


