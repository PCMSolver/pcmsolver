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
#include <stdlib.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "sas.h"
#include "eval.h"
#include "smooth.h"



void qacl_plas_wufg(double *t,int nb,int k,ns_curv1D *B)
{int i;
double *x;
x=(double *)malloc((nb+1)*sizeof(double));
for(i=0;i<=nb;i++)
	x[i]=0.0;
x[k]=1.0;
fepn_inte_fohk(t,x,nb,B);
free(x);
}



void mesv_inte_rujm(ns_surf S,int nb,int type,
double val_fix,ns_curv *F)
{int i;
double *t,step;
point *val;
t  =(double *)malloc((nb+1)*sizeof(double));
val=(point *)malloc((nb+1)*sizeof(point));
step=1.0/(double)nb;
for(i=0;i<=nb;i++)
	{t[i]=(double)i*step;
	if(type==1)
		cilj_eval_qelf(S,t[i],val_fix,&val[i]);
	if(type==2)
		cilj_eval_qelf(S,val_fix,t[i],&val[i]);
	}
ralc_inte_zuts(t,val,nb,F);
free(val);
free(t);
}



void gozn_prop_kupl(ns_surf S,int N,int M,ns_curv1D *phi,
ns_curv1D *psi,ns_curv *f,ns_curv *g,point **x)
{int i,j,k,l;
double *u,*v,step_u,step_v;

u=(double *)malloc((M+1)*sizeof(double));
v=(double *)malloc((N+1)*sizeof(double));
step_u=1.0/(double)M;
step_v=1.0/(double)N;
for(k=0;k<=M;k++)
	u[k]=(double)k*step_u;
for(l=0;l<=N;l++)
	v[l]=(double)l*step_v;
for(i=0;i<=M;i++)
	qacl_plas_wufg(u,M,i,&phi[i]);
for(j=0;j<=N;j++)
	qacl_plas_wufg(v,N,j,&psi[j]);

for(k=0;k<=M;k++)
for(l=0;l<=N;l++)
	cilj_eval_qelf(S,u[k],v[l],&x[k][l]);

for(j=0;j<=N;j++)
	mesv_inte_rujm(S,M,1,v[j],&f[j]);
for(i=0;i<=M;i++)
	mesv_inte_rujm(S,N,2,u[i],&g[i]);

free(u);
free(v);
}


void mosn_eval_gonr(double u,double v,int N,int M,
ns_curv1D *phi,ns_curv1D *psi,ns_curv *f,
ns_curv *g,point **x,point *sol)
{int i,j;
double coeff,*psi_val,*phi_val;
point term1,term2,term3,temp;

phi_val=(double *)malloc((M+1)*sizeof(double));
psi_val=(double *)malloc((N+1)*sizeof(double));
for(i=0;i<=M;i++)
	phi_val[i]=jurw_eval_voch(phi[i],u);
for(j=0;j<=N;j++)
	psi_val[j]=jurw_eval_voch(psi[j],v);

term1.absi=0.0;	term1.ordo=0.0;	term1.cote=0.0;
for(i=0;i<=M;i++)
	{cuwd_eval_nivk(g[i],v,&temp);
	coeff=phi_val[i];
	term1.absi=term1.absi+coeff*temp.absi;
	term1.ordo=term1.ordo+coeff*temp.ordo;
	term1.cote=term1.cote+coeff*temp.cote;
	}

term2.absi=0.0;	term2.ordo=0.0;	term2.cote=0.0;
for(j=0;j<=N;j++)
	{cuwd_eval_nivk(f[j],u,&temp);
	coeff=psi_val[j];
	term2.absi=term2.absi+coeff*temp.absi;
	term2.ordo=term2.ordo+coeff*temp.ordo;
	term2.cote=term2.cote+coeff*temp.cote;
	}

term3.absi=0.0;	term3.ordo=0.0;	term3.cote=0.0;
for(i=0;i<=M;i++)
for(j=0;j<=N;j++)
	{term3.absi=term3.absi+x[i][j].absi*phi_val[i]*psi_val[j];
	term3.ordo =term3.ordo+x[i][j].ordo*phi_val[i]*psi_val[j];
	term3.cote =term3.cote+x[i][j].cote*phi_val[i]*psi_val[j];
	}
free(phi_val);
free(psi_val);

sol->absi=term1.absi+term2.absi-term3.absi;
sol->ordo=term1.ordo+term2.ordo-term3.ordo;
sol->cote=term1.cote+term2.cote-term3.cote;
}


void duqn_extr_qepw(int m,ns_surf S,ns_curv *C)
{int i,k,n;
k=S.ku;
n=S.nu;
C->k=k;
C->n=n;
for(i=0;i<=n;i++)
	{getf_find_rogc_todj(S.d[i][m],&C->d[i]);
	C->w[i]=S.w[i][m];
	}
for(i=0;i<=n+k;i++)
	C->tau[i]=S.frknt[i];
C->v0=S.u0;
C->v1=S.u1;
C->prop3=1;
}


void capd_extr_widl(int m,ns_surf S,ns_curv *C)
{int i,j,k,n;
k=S.kv;
n=S.nv;
C->k=k;
C->n=n;
for(j=0;j<=n;j++)
	{getf_find_rogc_todj(S.d[m][j],&C->d[j]);
	C->w[j]=S.w[m][j];
	}
for(i=0;i<=n+k;i++)
	C->tau[i]=S.scknt[i];
C->v0=S.v0;
C->v1=S.v1;
C->prop3=1;
}


void tegj_smoo_gekc(int L_bic,int M_bic,ns_surf S_in,
int N_gord,int M_gord,ns_surf *S_out)
{int i,j;
double *u_bic,*v_bic,step_u,step_v;
ns_curv1D *phi,*psi;
ns_curv *f,*g,C_1,C_2,C_3,C_4;
prop_n_curv pnc_hrz,pnc_vrt;
prop_n_curv pnc_f,pnc_g,pnc_phi,pnc_psi;
point **x,**P;

pnc_f.k=4;	pnc_f.n=M_gord+2;
f=(ns_curv *)malloc((N_gord+1)*sizeof(ns_curv));
for(j=0;j<=N_gord;j++)
	foks_allo_vukp(pnc_f,&f[j]);
pnc_g.k=4;	pnc_g.n=N_gord+2;
g=(ns_curv *)malloc((M_gord+1)*sizeof(ns_curv));
for(i=0;i<=M_gord;i++)
	foks_allo_vukp(pnc_g,&g[i]);
x=(point **)malloc((M_gord+1)*sizeof(point));
for(i=0;i<=M_gord;i++)
	x[i]=(point *)malloc((N_gord+1)*sizeof(point));
pnc_phi.k=4;	pnc_phi.n=M_gord+2;
pnc_psi.k=4;	pnc_psi.n=N_gord+2;
phi=(ns_curv1D *)malloc((M_gord+1)*sizeof(ns_curv1D));
for(i=0;i<=M_gord;i++)
	fojb_allo_qugd(pnc_phi,&phi[i]);
psi=(ns_curv1D *)malloc((N_gord+1)*sizeof(ns_curv1D));
for(j=0;j<=N_gord;j++)
	fojb_allo_qugd(pnc_psi,&psi[j]);
gozn_prop_kupl(S_in,N_gord,M_gord,phi,psi,f,g,x);


u_bic=(double *)malloc((L_bic+1)*sizeof(double));
v_bic=(double *)malloc((M_bic+1)*sizeof(double));
step_u=1.0/(double)L_bic;
step_v=1.0/(double)M_bic;
for(i=0;i<=L_bic;i++)
	u_bic[i]=step_u*(double)i;
for(j=0;j<=M_bic;j++)
	v_bic[j]=step_v*(double)j;
P=(point **)malloc((L_bic+1)*sizeof(point));
for(i=0;i<=L_bic;i++)
	{P[i]=(point *)malloc((M_bic+1)*sizeof(point));
	for(j=0;j<=M_bic;j++)
		mosn_eval_gonr(u_bic[i],v_bic[j],N_gord,M_gord,phi,psi,f,g,x,&P[i][j]);
	}


pnc_hrz.k=S_in.ku;	pnc_hrz.n=S_in.nu;
pnc_vrt.k=S_in.kv;	pnc_vrt.n=S_in.nv;
foks_allo_vukp(pnc_hrz,&C_1);
foks_allo_vukp(pnc_hrz,&C_2);
foks_allo_vukp(pnc_vrt,&C_3);
foks_allo_vukp(pnc_vrt,&C_4);
duqn_extr_qepw(0,S_in,&C_1);
duqn_extr_qepw(S_in.nv,S_in,&C_2);
capd_extr_widl(0,S_in,&C_3);
capd_extr_widl(S_in.nu,S_in,&C_4);
reck_bicu_bavh(u_bic,v_bic,P,
L_bic,M_bic,C_1,C_2,C_3,C_4,S_out);
newt_dest_lefq(pnc_hrz,&C_1);
newt_dest_lefq(pnc_hrz,&C_2);
newt_dest_lefq(pnc_vrt,&C_3);
newt_dest_lefq(pnc_vrt,&C_4);

for(i=0;i<=M_gord;i++)
	free(x[i]);
free(x);
for(j=0;j<=N_gord;j++)
	newt_dest_lefq(pnc_f,&f[j]);
free(f);
for(i=0;i<=M_gord;i++)
	newt_dest_lefq(pnc_g,&g[i]);
free(g);
for(i=0;i<=M_gord;i++)
	gojt_dest_wujl(&phi[i]);
free(phi);	
for(j=0;j<=N_gord;j++)
	gojt_dest_wujl(&psi[j]);
free(psi);

for(i=0;i<=L_bic;i++)
	free(P[i]);
free(P);
free(u_bic);
free(v_bic);
}

 


void rusw_smoo_cohs(int N_gord,int M_gord,ns_surf *S)
{int L_bic,M_bic;
ns_surf S_temp;
prop_n_surf pns;
L_bic=S->nu-2;	
M_bic=S->nv-2;
pns.ku=4;		pns.kv=4;
pns.nu=S->nu;	pns.nv=S->nv;
juvm_allo_sehv(pns,&S_temp);
tegj_smoo_gekc(L_bic,M_bic,*S,N_gord,M_gord,&S_temp);
qucp_find_pogc_gecz(S_temp,S);
destroy_nurbs_surface_alx(pns,&S_temp);
}


void wemh_norm_mocl(ns_surf S,
double u,double v,double h,vect3D *N)
{point val,val_ph;
vect3D U,V;
cilj_eval_qelf(S,u,v,&val);

cilj_eval_qelf(S,u+h,v,&val_ph);
U.absi=(val_ph.absi-val.absi)/h;
U.ordo=(val_ph.ordo-val.ordo)/h;
U.cote=(val_ph.cote-val.cote)/h;

cilj_eval_qelf(S,u,v+h,&val_ph);
V.absi=(val_ph.absi-val.absi)/h;
V.ordo=(val_ph.ordo-val.ordo)/h;
V.cote=(val_ph.cote-val.cote)/h;

cofz_cros_fits(U,V,N);
qubr_norm_foqk(N);
}


void gurn_norm_tegm(ns_surf S,double u,
double v,double h,vect3D *N)
{double u_aux,v_aux,eps=1.0e-3;
if(u+h>1.0)
	u_aux=u-(1.0+eps)*h;
else
	u_aux=u;

if(v+h>1.0)
	v_aux=v-(1.0+eps)*h;
else
	v_aux=v;
wemh_norm_mocl(S,u_aux,v_aux,h,N);
}



void qumc_cent_hivf(ns_surf S,int side,ns_curv *C)
{int i,n;
double *t,step,u,v;
point *X,val,val_ph;
double h=1.0e-2;
n=S.nu-2;
t=(double *)malloc((n+1)*sizeof(double));
X=(point *)malloc((n+1)*sizeof(point));
step=1.0/(double)n;
for(i=0;i<=n;i++)
	{t[i]=(double)i*step;
	if(side==0)
		{u=t[i];
		cilj_eval_qelf(S,u,0.0,&val);
		cilj_eval_qelf(S,u,h,&val_ph);
		X[i].absi=(val_ph.absi-val.absi)/h;
		X[i].ordo=(val_ph.ordo-val.ordo)/h;
		X[i].cote=(val_ph.cote-val.cote)/h;
		}
	if(side==1)
		{v=t[i];
		cilj_eval_qelf(S,1.0,v,&val);
		cilj_eval_qelf(S,1.0-h,v,&val_ph);
		X[i].absi=(val_ph.absi-val.absi)/h;
		X[i].ordo=(val_ph.ordo-val.ordo)/h;
		X[i].cote=(val_ph.cote-val.cote)/h;
		}
	if(side==2)
		{u=t[i];
		cilj_eval_qelf(S,u,1.0,&val);
		cilj_eval_qelf(S,u,1.0-h,&val_ph);
		X[i].absi=(val_ph.absi-val.absi)/h;
		X[i].ordo=(val_ph.ordo-val.ordo)/h;
		X[i].cote=(val_ph.cote-val.cote)/h;
		}
	if(side==3)
		{v=t[i];
		cilj_eval_qelf(S,0.0,v,&val);
		cilj_eval_qelf(S,h,v,&val_ph);
		X[i].absi=(val_ph.absi-val.absi)/h;
		X[i].ordo=(val_ph.ordo-val.ordo)/h;
		X[i].cote=(val_ph.cote-val.cote)/h;
		}
	}
ralc_inte_zuts(t,X,n,C);
free(t);
free(X);
}



void gidc_boun_lirm(ns_surf S,int side,ns_curv *C)
{int i,n;
double *t,step,u,v;
point *X;
double h=1.0e-2;
n=S.nu-2;
t=(double *)malloc((n+1)*sizeof(double));
X=(point *)malloc((n+1)*sizeof(point));
step=1.0/(double)n;
for(i=0;i<=n;i++)
	{t[i]=(double)i*step;
	if(side==0)
		{u=t[i];
		gurn_norm_tegm(S,u,0.0,h,&X[i]);
		}
	if(side==1)
		{v=t[i];
		gurn_norm_tegm(S,1.0,v,h,&X[i]);
		}
	if(side==2)
		{u=t[i];
		gurn_norm_tegm(S,u,1.0,h,&X[i]);
		}
	if(side==3)
		{v=t[i];
		gurn_norm_tegm(S,0.0,v,h,&X[i]);
		}
	}
ralc_inte_zuts(t,X,n,C);
free(t);
free(X);
}


double qofs_find_gehk_nurk(double delta,double t)
{double res,*b,*c,t_tilde;
if((delta>=1.0)||(delta<=0.0))
	{fprintf(tmpout,"Awry value of delta\n");
	exit(0);
	}
if(t<=delta)
	{b=(double *)malloc(4*sizeof(double));
	b[0]=0.0;	b[1]=0.0;
	b[2]=1.0;	b[3]=1.0;
	t_tilde=delta*t;
	res=mijt_scal_letp(4,b,t_tilde);
	free(b);
	}
else if(t>=1.0-delta)
	{c=(double *)malloc(4*sizeof(double));
	c[0]=1.0;	c[1]=1.0;
	c[2]=0.0;	c[3]=0.0;
	t_tilde=(t+delta-1.0)/delta;
	res=mijt_scal_letp(4,c,t_tilde);
	free(c);
	}
else
	res=1.0;
return res;
}



double humr_blen_rols(double u,double v)
{double lambda,delta=0.1,fact1,fact2;
fact1=qofs_find_gehk_nurk(delta,u);
fact2=qofs_find_gehk_nurk(delta,v);
lambda=fact1*fact2;
return lambda;
}


void vijb_take_quwf(point temp1,point temp2,point *X)
{X->absi=0.5*(temp1.absi+temp2.absi);
X->ordo=0.5*(temp1.ordo+temp2.ordo);
X->cote=0.5*(temp1.cote+temp2.cote);
}



void cuvk_find_duhr_desl(point **P,int L,int M,ns_curv A,ns_curv B,
ns_curv C,ns_curv D)
{int i,j;
double dist,step_u,step_v,u,v;
vect3D dir;
point temp1,temp2;

step_u=1.0/(double)L;
for(i=2;i<=L-2;i++)
	{dist=wodt_dist_gilq(P[i][0],P[i][1]);
	u=(double)i*step_u;
	cuwd_eval_nivk(A,u,&dir);
	qubr_norm_foqk(&dir);
	P[i][1].absi=P[i][0].absi+dist*dir.absi;
	P[i][1].ordo=P[i][0].ordo+dist*dir.ordo;
	P[i][1].cote=P[i][0].cote+dist*dir.cote;
	}

for(i=2;i<=L-2;i++)
	{dist=wodt_dist_gilq(P[i][M],P[i][M-1]);
	u=(double)i*step_u;
	cuwd_eval_nivk(C,u,&dir);
	qubr_norm_foqk(&dir);
	P[i][M-1].absi=P[i][M].absi-dist*dir.absi;
	P[i][M-1].ordo=P[i][M].ordo-dist*dir.ordo;
	P[i][M-1].cote=P[i][M].cote-dist*dir.cote;
	}

step_v=1.0/(double)M;
for(j=2;j<=M-2;j++)
	{dist=wodt_dist_gilq(P[0][j],P[1][j]);
	v=(double)j*step_v;
	cuwd_eval_nivk(D,v,&dir);
	qubr_norm_foqk(&dir);
	P[1][j].absi=P[0][j].absi+dist*dir.absi;
	P[1][j].ordo=P[0][j].ordo+dist*dir.ordo;
	P[1][j].cote=P[0][j].cote+dist*dir.cote;
	}

for(j=2;j<=M-2;j++)
	{dist=wodt_dist_gilq(P[L][j],P[L-1][j]);
	v=(double)j*step_v;
	cuwd_eval_nivk(B,v,&dir);
	qubr_norm_foqk(&dir);
	P[L-1][j].absi=P[L][j].absi-dist*dir.absi;
	P[L-1][j].ordo=P[L][j].ordo-dist*dir.ordo;
	P[L-1][j].cote=P[L][j].cote-dist*dir.cote;
	}

dist=wodt_dist_gilq(P[1][0],P[1][1]);
u=step_u;
cuwd_eval_nivk(A,u,&dir);
qubr_norm_foqk(&dir);
temp1.absi=P[1][0].absi+dist*dir.absi;
temp1.ordo=P[1][0].ordo+dist*dir.ordo;
temp1.cote=P[1][0].cote+dist*dir.cote;

dist=wodt_dist_gilq(P[0][1],P[1][1]);
v=step_v;
cuwd_eval_nivk(D,v,&dir);
qubr_norm_foqk(&dir);
temp2.absi=P[0][1].absi+dist*dir.absi;
temp2.ordo=P[0][1].ordo+dist*dir.ordo;
temp2.cote=P[0][1].cote+dist*dir.cote;

vijb_take_quwf(temp1,temp2,&P[1][1]);

dist=wodt_dist_gilq(P[L-1][0],P[L-1][1]);
u=((double)L-1.0)*step_u;
cuwd_eval_nivk(A,u,&dir);
qubr_norm_foqk(&dir);
temp1.absi=P[L-1][0].absi+dist*dir.absi;
temp1.ordo=P[L-1][0].ordo+dist*dir.ordo;
temp1.cote=P[L-1][0].cote+dist*dir.cote;

dist=wodt_dist_gilq(P[L][1],P[L-1][1]);
v=step_v;
cuwd_eval_nivk(B,v,&dir);
qubr_norm_foqk(&dir);
temp2.absi=P[L][1].absi-dist*dir.absi;
temp2.ordo=P[L][1].ordo-dist*dir.ordo;
temp2.cote=P[L][1].cote-dist*dir.cote;

vijb_take_quwf(temp1,temp2,&P[L-1][1]);

dist=wodt_dist_gilq(P[1][M],P[1][M-1]);
u=step_u;
cuwd_eval_nivk(C,u,&dir);
qubr_norm_foqk(&dir);
temp1.absi=P[1][M].absi-dist*dir.absi;
temp1.ordo=P[1][M].ordo-dist*dir.ordo;
temp1.cote=P[1][M].cote-dist*dir.cote;

dist=wodt_dist_gilq(P[0][M-1],P[1][M-1]);
v=((double)M-1.0)*step_v;
cuwd_eval_nivk(D,v,&dir);
qubr_norm_foqk(&dir);
temp2.absi=P[0][M-1].absi+dist*dir.absi;
temp2.ordo=P[0][M-1].ordo+dist*dir.ordo;
temp2.cote=P[0][M-1].cote+dist*dir.cote;

vijb_take_quwf(temp1,temp2,&P[1][M-1]);

dist=wodt_dist_gilq(P[L-1][M],P[L-1][M-1]);
u=((double)L-1.0)*step_u;
cuwd_eval_nivk(C,u,&dir);
qubr_norm_foqk(&dir);
temp1.absi=P[L-1][M].absi-dist*dir.absi;
temp1.ordo=P[L-1][M].ordo-dist*dir.ordo;
temp1.cote=P[L-1][M].cote-dist*dir.cote;

dist=wodt_dist_gilq(P[L][M-1],P[L-1][M-1]);
v=((double)M-1.0)*step_v;
cuwd_eval_nivk(B,v,&dir);
qubr_norm_foqk(&dir);
temp2.absi=P[L][M-1].absi-dist*dir.absi;
temp2.ordo=P[L][M-1].ordo-dist*dir.ordo;
temp2.cote=P[L][M-1].cote-dist*dir.cote;

vijb_take_quwf(temp1,temp2,&P[L-1][M-1]);
}


void hocr_boun_hecr(int L_bic,int M_bic,ns_surf S_in,
ns_curv A,ns_curv B,ns_curv C,ns_curv D,
ns_surf *S_out)
{int i,j;
double *u_bic,*v_bic,step_u,step_v;
point **P;
ns_curv C_1,C_2,C_3,C_4;
prop_n_curv pnc_hrz,pnc_vrt;
u_bic=(double *)malloc((L_bic+1)*sizeof(double));
v_bic=(double *)malloc((M_bic+1)*sizeof(double));
step_u=1.0/(double)L_bic;
step_v=1.0/(double)M_bic;
for(i=0;i<=L_bic;i++)
	u_bic[i]=step_u*(double)i;
for(j=0;j<=M_bic;j++)
	v_bic[j]=step_v*(double)j;
P=(point **)malloc((L_bic+1)*sizeof(point));
for(i=0;i<=L_bic;i++)
	{P[i]=(point *)malloc((M_bic+1)*sizeof(point));
	for(j=0;j<=M_bic;j++)
		cilj_eval_qelf(S_in,u_bic[i],v_bic[j],&P[i][j]);
	}
cuvk_find_duhr_desl(P,L_bic,M_bic,A,B,C,D);

pnc_hrz.k=S_in.ku;	pnc_hrz.n=S_in.nu;
pnc_vrt.k=S_in.kv;	pnc_vrt.n=S_in.nv;
foks_allo_vukp(pnc_hrz,&C_1);
foks_allo_vukp(pnc_hrz,&C_2);
foks_allo_vukp(pnc_vrt,&C_3);
foks_allo_vukp(pnc_vrt,&C_4);
duqn_extr_qepw(0,S_in,&C_1);
duqn_extr_qepw(S_in.nv,S_in,&C_2);
capd_extr_widl(0,S_in,&C_3);
capd_extr_widl(S_in.nu,S_in,&C_4);
reck_bicu_bavh(u_bic,v_bic,P,L_bic,
M_bic,C_1,C_2,C_3,C_4,S_out);
newt_dest_lefq(pnc_hrz,&C_1);
newt_dest_lefq(pnc_hrz,&C_2);
newt_dest_lefq(pnc_vrt,&C_3);
newt_dest_lefq(pnc_vrt,&C_4);

for(i=0;i<=L_bic;i++)
	free(P[i]);
free(P);
free(u_bic);
free(v_bic);
}



void gucd_mixt_vemg(int L_bic,int M_bic,ns_surf S_in_1,
ns_surf S_in_2,ns_surf *S_out)
{int i,j;
double *u_bic,*v_bic,step_u,step_v,lambda;
point **P,temp1,temp2;
ns_curv C_1,C_2,C_3,C_4;
prop_n_curv pnc_hrz,pnc_vrt;
u_bic=(double *)malloc((L_bic+1)*sizeof(double));
v_bic=(double *)malloc((M_bic+1)*sizeof(double));
step_u=1.0/(double)L_bic;
step_v=1.0/(double)M_bic;
for(i=0;i<=L_bic;i++)
	u_bic[i]=step_u*(double)i;
for(j=0;j<=M_bic;j++)
	v_bic[j]=step_v*(double)j;
P=(point **)malloc((L_bic+1)*sizeof(point));
for(i=0;i<=L_bic;i++)
	{P[i]=(point *)malloc((M_bic+1)*sizeof(point));
	for(j=0;j<=M_bic;j++)
		{cilj_eval_qelf(S_in_1,u_bic[i],v_bic[j],&temp1);
		cilj_eval_qelf(S_in_2,u_bic[i],v_bic[j],&temp2);
		lambda=humr_blen_rols(u_bic[i],v_bic[j]);
		P[i][j].absi=lambda*temp2.absi+(1.0-lambda)*temp1.absi;
		P[i][j].ordo=lambda*temp2.ordo+(1.0-lambda)*temp1.ordo;
		P[i][j].cote=lambda*temp2.cote+(1.0-lambda)*temp1.cote;
		}
	}

pnc_hrz.k=S_in_1.ku;	pnc_hrz.n=S_in_1.nu;
pnc_vrt.k=S_in_1.kv;	pnc_vrt.n=S_in_1.nv;
foks_allo_vukp(pnc_hrz,&C_1);  
foks_allo_vukp(pnc_hrz,&C_2);
foks_allo_vukp(pnc_vrt,&C_3);
foks_allo_vukp(pnc_vrt,&C_4);
duqn_extr_qepw(0,S_in_1,&C_1);
duqn_extr_qepw(S_in_1.nv,S_in_1,&C_2);
capd_extr_widl(0,S_in_1,&C_3);
capd_extr_widl(S_in_1.nu,S_in_1,&C_4);
reck_bicu_bavh(u_bic,v_bic,P,L_bic,
M_bic,C_1,C_2,C_3,C_4,S_out);
newt_dest_lefq(pnc_hrz,&C_1);
newt_dest_lefq(pnc_hrz,&C_2);
newt_dest_lefq(pnc_vrt,&C_3);
newt_dest_lefq(pnc_vrt,&C_4);

for(i=0;i<=L_bic;i++)
	free(P[i]);
free(P);
free(u_bic);
free(v_bic);
}


void mups_corn_cogj(int side,ns_surf S,
point *A,point *B)
{int n_u,n_v;
n_u=S.nu;
n_v=S.nv;
switch(side)
	{case 0:
		getf_find_rogc_todj(S.d[0][0],A);
		getf_find_rogc_todj(S.d[n_u][0],B);
		break;
	case 1:
		getf_find_rogc_todj(S.d[n_u][0],A);
		getf_find_rogc_todj(S.d[n_u][n_v],B);
		break;
	case 2:
		getf_find_rogc_todj(S.d[n_u][n_v],A);
		getf_find_rogc_todj(S.d[0][n_v],B);
		break;
	case 3:
		getf_find_rogc_todj(S.d[0][0],A);
		getf_find_rogc_todj(S.d[0][n_v],B);
		break;
	}
}



void humz_find_semv_nuwc(vect3D W_cur,vect3D W_ngb,
vect3D N_cur,vect3D N_ngb,vect3D *sol)
{double c1,c2;
vect3D N,E1,E2;

N.absi=0.5*(N_cur.absi+N_ngb.absi);
N.ordo=0.5*(N_cur.ordo+N_ngb.ordo);
N.cote=0.5*(N_cur.cote+N_ngb.cote);
qubr_norm_foqk(&N);
cofz_cros_fits(N,W_cur,&E1);
qubr_norm_foqk(&E1);
cofz_cros_fits(N,E1,&E2);
qubr_norm_foqk(&E2);

c1=rocv_scal_toqc(E1,W_cur);
c2=rocv_scal_toqc(E2,W_cur);
sol->absi=c1*E1.absi+c2*E2.absi;
sol->ordo=c1*E1.ordo+c2*E2.ordo;
sol->cote=c1*E1.cote+c2*E2.cote;
qubr_norm_foqk(sol);
}



void wuzv_esti_weql(ns_surf *S,int id_cur,int side_cur,
int id_ngb,ns_curv **T,ns_curv **N,ns_curv *t)
{int side_ngb,i;
double sml,d1,d2,dis;
point A_cur,B_cur,A_ngb,B_ngb;
prop_n_curv pnc;
ns_curv t_cur,t_ngb,n_cur,n_ngb;
vect3D temp;

pnc.k=S[id_cur].ku;
pnc.n=S[id_cur].nu;
foks_allo_vukp(pnc,&t_cur);
foks_allo_vukp(pnc,&t_ngb);
foks_allo_vukp(pnc,&n_cur);
foks_allo_vukp(pnc,&n_ngb);
mups_corn_cogj(side_cur,S[id_cur],&A_cur,&B_cur);
sml=LARGE_NUMBER;
for(i=0;i<4;i++)
	{mups_corn_cogj(i,S[id_ngb],&A_ngb,&B_ngb);
	d1=wodt_dist_gilq(A_cur,A_ngb)+wodt_dist_gilq(B_cur,B_ngb);
	d2=wodt_dist_gilq(B_cur,A_ngb)+wodt_dist_gilq(A_cur,B_ngb);
	if(d1<d2)	dis=d1;
	else		dis=d2;	
	if(dis<sml)
		{sml=dis;
		side_ngb=i;
		}
	}


if((side_cur==0)||(side_cur==3))
	{zobm_find_wumq_kihf(T[id_cur][side_cur],&t_cur);
	zobm_find_wumq_kihf(T[id_ngb][side_ngb],&t_ngb);
	for(i=0;i<=t_ngb.n;i++)
		{t_ngb.d[i].absi=-t_ngb.d[i].absi;
		t_ngb.d[i].ordo =-t_ngb.d[i].ordo;
		t_ngb.d[i].cote =-t_ngb.d[i].cote;
		}
	}
if((side_cur==1)||(side_cur==2))
	{zobm_find_wumq_kihf(T[id_cur][side_cur],&t_cur);
	for(i=0;i<=t_cur.n;i++)
		{t_cur.d[i].absi=-t_cur.d[i].absi;
		t_cur.d[i].ordo =-t_cur.d[i].ordo;
		t_cur.d[i].cote =-t_cur.d[i].cote;
		}
	zobm_find_wumq_kihf(T[id_ngb][side_ngb],&t_ngb);
	}
pobd_flip_kejt(&t_ngb);

zobm_find_wumq_kihf(N[id_cur][side_cur],&n_cur);
zobm_find_wumq_kihf(N[id_ngb][side_ngb],&n_ngb);
zobm_find_wumq_kihf(t_cur,t);
for(i=0;i<=t_cur.n;i++)
	{humz_find_semv_nuwc(t_cur.d[i],t_ngb.d[i],n_cur.d[i],n_ngb.d[i],&temp);
	getf_find_rogc_todj(temp,&t->d[i]);	
	}
newt_dest_lefq(pnc,&t_cur);
newt_dest_lefq(pnc,&t_ngb);
newt_dest_lefq(pnc,&n_cur);
newt_dest_lefq(pnc,&n_ngb);
}


double fans_dist_qefr(point A,point B,ns_surf S)
{int i;
double d1,d2,dis,res;
point a,b;
res=LARGE_NUMBER;
for(i=0;i<4;i++)
	{mups_corn_cogj(i,S,&a,&b);
	d1=wodt_dist_gilq(a,A)+wodt_dist_gilq(b,B);
	d2=wodt_dist_gilq(b,A)+wodt_dist_gilq(a,B);
	if(d1<d2)	dis=d1;
	else		dis=d2;
	if(dis<res)
		res=dis;
	}
return res;
}



void legs_side_dolp(ns_surf *S,fajor_sion3D QUAD,
int id_cur,int side_cur,ns_curv **T,ns_curv **N,ns_curv *t)
{int id_ngb,ed[4],i,E,e1,e2,f;
double dis,sml;
point A,B;
mups_corn_cogj(side_cur,S[id_cur],&A,&B);
ed[0]=QUAD.elem[id_cur].frkt;
ed[1]=QUAD.elem[id_cur].sckt;
ed[2]=QUAD.elem[id_cur].trkt;
ed[3]=QUAD.elem[id_cur].ftkt;
sml=LARGE_NUMBER;
for(i=0;i<4;i++)
	{E=ed[i];
	e1=QUAD.kt[E].frent;
	e2=QUAD.kt[E].scent;
	if(e1==id_cur)	f=e2;
	if(e2==id_cur)	f=e1;
	dis=fans_dist_qefr(A,B,S[f]);
	if(dis<sml)
		{sml=dis;
		id_ngb=f;
		}
	}
wuzv_esti_weql(S,id_cur,side_cur,id_ngb,T,N,t);
}


void tofd_find_sofq_fovr(ns_surf *S,fajor_sion3D QUAD,
int id_cur,ns_curv **T,ns_curv **N)
{int L_bic,M_bic;
prop_n_surf pns;
prop_n_curv pnc_hrz,pnc_vrt;
ns_curv A,B,C,D;
ns_surf S_1,S_2;

pns.ku=S[id_cur].ku;
pns.kv=S[id_cur].kv;
pns.nu=S[id_cur].nu;
pns.nv=S[id_cur].nv;
juvm_allo_sehv(pns,&S_1);
juvm_allo_sehv(pns,&S_2);
pnc_hrz.k=S[id_cur].ku;		pnc_hrz.n=S[id_cur].nu;
pnc_vrt.k=S[id_cur].kv;		pnc_vrt.n=S[id_cur].nv;

foks_allo_vukp(pnc_hrz,&A);
foks_allo_vukp(pnc_vrt,&B);
foks_allo_vukp(pnc_hrz,&C);
foks_allo_vukp(pnc_vrt,&D);
legs_side_dolp(S,QUAD,id_cur,0,T,N,&A);
legs_side_dolp(S,QUAD,id_cur,1,T,N,&B);
legs_side_dolp(S,QUAD,id_cur,2,T,N,&C);
legs_side_dolp(S,QUAD,id_cur,3,T,N,&D);

L_bic=S[id_cur].nu-2;	
M_bic=S[id_cur].nv-2;
hocr_boun_hecr(L_bic,M_bic,S[id_cur],A,B,C,D,&S_1);
qucp_find_pogc_gecz(S[id_cur],&S_2);
gucd_mixt_vemg(L_bic,M_bic,S_1,S_2,&S[id_cur]);
destroy_nurbs_surface_alx(pns,&S_2);
destroy_nurbs_surface_alx(pns,&S_1);
} 


void huvg_find_mucv_newr(fajor_sion3D QUAD,
ns_surf *surf,int n_surf)
{int i,j;
ns_curv **T,**N;
prop_n_curv pnc;

pnc.k=4;
pnc.n=surf[0].nu;
T=(ns_curv **)malloc(n_surf*sizeof(ns_curv));
N=(ns_curv **)malloc(n_surf*sizeof(ns_curv));
for(i=0;i<n_surf;i++)
	{T[i]=(ns_curv *)malloc(4*sizeof(ns_curv));
	N[i]=(ns_curv *)malloc(4*sizeof(ns_curv));
	for(j=0;j<4;j++)
		{foks_allo_vukp(pnc,&T[i][j]);
		foks_allo_vukp(pnc,&N[i][j]);
		}
	}
for(i=0;i<n_surf;i++)
for(j=0;j<4;j++)
	{qumc_cent_hivf(surf[i],j,&T[i][j]);
	gidc_boun_lirm(surf[i],j,&N[i][j]);
	}

for(i=0;i<n_surf;i++)
	{tofd_find_sofq_fovr(surf,QUAD,i,T,N);
	
	}

for(i=0;i<n_surf;i++)
	{for(j=0;j<4;j++)
		{newt_dest_lefq(pnc,&T[i][j]);
		newt_dest_lefq(pnc,&N[i][j]);
		}
	free(T[i]);
	free(N[i]);
	}
free(T);
free(N);
}
 
