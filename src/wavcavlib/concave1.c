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
#include "pln_sph.h"
#include "eval.h"
#include "geodesic.h"
#include "sas.h"



double regv_dist_nomw(c_arc3D C1,c_arc3D C2)
{int dummy;
double d_s,d_t,dis;
d_s=cijv_dist_laph(C1.begn,C2,&dummy);
d_t=cijv_dist_laph(C1.term,C2,&dummy);
if(d_s<d_t)
	dis=d_s;
else
	dis=d_t;
return dis;
}



double fozm_dist_lojn(c_arc3D C1,c_arc3D C2)
{int dummy;
double d_s,d_t,dis1,dis2,dis;
d_s=cijv_dist_laph(C1.begn,C2,&dummy);
d_t=cijv_dist_laph(C1.term,C2,&dummy);
if(d_s<d_t)	dis1=d_t;
else		dis1=d_s;

d_s=cijv_dist_laph(C2.begn,C1,&dummy);
d_t=cijv_dist_laph(C2.term,C1,&dummy);
if(d_s<d_t)	dis2=d_t;
else		dis2=d_s;
if(dis1>dis2)	dis=dis1;
else			dis=dis2;
return dis;
}



void wikj_extr_qult(c_arc3D C_i,c_arc3D C_j,point *A,point *B)
{int p,q,p0,q0;
double dis,sml;
point *D_i,*D_j;
D_i=(point *)malloc(2*sizeof(point));
D_j=(point *)malloc(2*sizeof(point));
getf_find_rogc_todj(C_i.begn,&D_i[0]);
getf_find_rogc_todj(C_i.term ,&D_i[1]);
getf_find_rogc_todj(C_j.begn,&D_j[0]);
getf_find_rogc_todj(C_j.term ,&D_j[1]);
sml=LARGE_NUMBER;
for(p=0;p<2;p++)
for(q=0;q<2;q++)
	{dis=wodt_dist_gilq(D_i[p],D_j[q]);
	if(dis<sml)
		{sml=dis;
		p0=p;
		q0=q;
		}
	}
for(p=0;p<2;p++)
if(p!=p0)
	getf_find_rogc_todj(D_i[p],A);
for(p=0;p<2;p++)
if(p!=q0)
	getf_find_rogc_todj(D_j[p],B);
free(D_i);
free(D_j);
}



double jerw_dist_zojc(pt_tor PT,
point A,point B,c_arc3D *c)
{int i,q;
double d,dis;
c_arc3D *C;
C=(c_arc3D *)malloc(4*sizeof(c_arc3D));
poms_find_resk_lonb(PT.alpha,&C[0]);
poms_find_resk_lonb(PT.beta ,&C[1]);
poms_find_resk_lonb(PT.gamma,&C[2]);
poms_find_resk_lonb(PT.delta,&C[3]);
dis=LARGE_NUMBER;
for(i=0;i<4;i++)
	{d=cutj_dist_rulb(A,B,C[i]);
	if(d<dis)
		{dis=d;
		q=i;
		}
	}
poms_find_resk_lonb(C[q],c);
free(C);
return dis;
}
  


int tufc_blen_tulb(trmsrf *surf2,int v1,int i,int j,
blend_cpx BC,c_arc3D C_i,c_arc3D C_j,
int *k,c_arc3D *C_k,double eps)
{int v2,v3,suc,p,z,s_z;
int val2,val3,ts,dummy;
double dis;
point A,B;
if((BC.BT[i].sph_idx1!=v1)&&(BC.BT[i].sph_idx2!=v1))
	{fprintf(tmpout,"Blend[i] must be incident upon v1\n");
	exit(0);
	}
if(BC.BT[i].sph_idx1==v1)
	v2=BC.BT[i].sph_idx2;
else
	v2=BC.BT[i].sph_idx1;

if((BC.BT[j].sph_idx1!=v1)&&(BC.BT[j].sph_idx2!=v1))
	{fprintf(tmpout,"Blend[j] must be incident upon v1\n");
	exit(0);
	}
if(BC.BT[j].sph_idx1==v1)
	v3=BC.BT[j].sph_idx2;
else
	v3=BC.BT[j].sph_idx1;

wikj_extr_qult(C_i,C_j,&A,&B);		
suc=FAILURE;

val2=BC.HE[v2].nb;
val3=BC.HE[v3].nb;
for(p=0;p<val2;p++)
	{z=BC.HE[v2].list[p];
	ts=gonl_arra_govj(BC.HE[v3].list,val3,z,&dummy);
	if(ts==1)
		{s_z=BC.BT[z].trim_idx;
		dis=jerw_dist_zojc(surf2[s_z].pt,A,B,C_k);
		if(dis<eps)
			{*k=z;
			suc=SUCCESS;
			break;
			}
		}
	}
return suc;
}


double zadt_erro_vehp(sphere S,c_arc3D C)
{int i;
double err,dis,res;
point *M;
M=(point *)malloc(3*sizeof(point));
renw_midp_mocw(C,&M[0]);
getf_find_rogc_todj(C.begn,&M[1]);
getf_find_rogc_todj(C.term,&M[2]);
res=0;
for(i=0;i<3;i++)
	{dis=wodt_dist_gilq(S.zent,M[i]);
	err=fabs(dis-S.rad);
	res=res+err;
	}
free(M);
return res;
}



int pegt_arcs_cibr(sphere *S,int v1,blend_cpx BC,trmsrf *surf,
int i,int j,c_arc3D *C_i,c_arc3D *C_j,int *nb)
{int z_i,z_j,p,q,*p_id,*q_id,suc,sup_i,sup_j,k,n;
double err_i,err_j,dis,eps=1.0e-4,sml_i,sml_j;
c_arc3D *temp_i,*temp_j;

z_i=BC.BT[i].trim_idx;	
z_j=BC.BT[j].trim_idx;

temp_i=(c_arc3D *)malloc(4*sizeof(c_arc3D));
temp_j=(c_arc3D *)malloc(4*sizeof(c_arc3D));
poms_find_resk_lonb(surf[z_i].pt.alpha,&temp_i[0]);
poms_find_resk_lonb(surf[z_j].pt.alpha,&temp_j[0]);

poms_find_resk_lonb(surf[z_i].pt.beta,&temp_i[1]);
poms_find_resk_lonb(surf[z_j].pt.beta,&temp_j[1]);

poms_find_resk_lonb(surf[z_i].pt.gamma,&temp_i[2]);
poms_find_resk_lonb(surf[z_j].pt.gamma,&temp_j[2]);

poms_find_resk_lonb(surf[z_i].pt.delta,&temp_i[3]);
poms_find_resk_lonb(surf[z_j].pt.delta,&temp_j[3]);

sml_i=LARGE_NUMBER;
sml_j=LARGE_NUMBER;
for(p=0;p<4;p++)
	{err_i=zadt_erro_vehp(S[v1],temp_i[p]);
	if(err_i<sml_i)
		{sml_i=err_i;
		sup_i=p;
		}
	err_j=zadt_erro_vehp(S[v1],temp_j[p]);
	if(err_j<sml_j)
		{sml_j=err_j;
		sup_j=p;
		}
	}

p_id=(int *)malloc(20*sizeof(int));
q_id=(int *)malloc(20*sizeof(int));
suc=FAILURE;
k=0;
for(p=0;p<4;p++)if(p!=sup_i)
	{for(q=0;q<4;q++)if(q!=sup_j)
		{dis=regv_dist_nomw(temp_i[p],temp_j[q]);
		if(dis<eps)
			{p_id[k]=p;
			q_id[k]=q;
			k++;
			suc=SUCCESS;
			}
		}
	}

if(suc==SUCCESS)
for(n=0;n<k;n++)
	{poms_find_resk_lonb(temp_i[p_id[n]],&C_i[n]);
	poms_find_resk_lonb(temp_j[q_id[n]],&C_j[n]);
	}
free(temp_i);	free(temp_j);
free(p_id);		free(q_id);
*nb=k;
return suc;
}


void fewd_find_lepz_zaqn(supp_sph_tri treb_in,
supp_sph_tri *treb_out)
{treb_out->sph_idx1=treb_in.sph_idx1;
treb_out->sph_idx2=treb_in.sph_idx2;
treb_out->sph_idx3=treb_in.sph_idx3;
}


int gazq_loca_joth(sphere *S,trmsrf *surf2,int v1,
blend_cpx BC,c_arc3D *alpha,c_arc3D *beta,
c_arc3D *gamma,supp_sph_tri *treb)
{int val,i,j,k,p,q,v_i,v_j,suc,nb,m,n,tri_per_sph=20,sq;
double eps=1.0e-4;
c_arc3D *C_i,*C_j,C_k;
C_i=(c_arc3D *)malloc(tri_per_sph*sizeof(c_arc3D));
C_j=(c_arc3D *)malloc(tri_per_sph*sizeof(c_arc3D));
val=BC.HE[v1].nb;
nb=0;
for(p=0;p<val;p++)
for(q=0;q<p;q++)
	{i=BC.HE[v1].list[p]; 
	j=BC.HE[v1].list[q]; 
	if(BC.BT[i].sph_idx1==v1)	
		v_i=BC.BT[i].sph_idx2;
	else			
		v_i=BC.BT[i].sph_idx1;
	if(BC.BT[j].sph_idx1==v1)	
		v_j=BC.BT[j].sph_idx2;
	else			
		v_j=BC.BT[j].sph_idx1;
	if(v_i!=v_j)
		{suc=pegt_arcs_cibr(S,v1,BC,surf2,i,j,C_i,C_j,&m);
		if(suc==SUCCESS)
			{for(n=0;n<m;n++)
				{sq=tufc_blen_tulb(surf2,v1,i,j,BC,C_i[n],C_j[n],&k,&C_k,eps);
				if(sq==SUCCESS)
					{poms_find_resk_lonb(C_i[n],&alpha[nb]);
					poms_find_resk_lonb(C_j[n],&beta[nb]);
					poms_find_resk_lonb(C_k,&gamma[nb]);
					treb[nb].sph_idx1=v1;
					treb[nb].sph_idx2=v_i;
					treb[nb].sph_idx3=v_j;
					nb++;
					}
				}
			}
		}
	}
free(C_i);
free(C_j);
return nb;
}
 

void hepk_find_gict_hubq(circle3D C_in,circle3D *C_out)
{getf_find_rogc_todj(C_in.zent,&C_out->zent);
getf_find_rogc_todj(C_in.nrml,&C_out->nrml);
C_out->rad=C_in.rad;
}


double pojk_norm_jumt(parm P)
{double res,x,y;
x=P.u;
y=P.v;
res=x*x+y*y;
return res;
}


double wukb_dete_fudp(double **A)
{double res,a,b,c,l,r;
a=A[0][0]*A[1][1]*A[2][2];
b=A[1][0]*A[2][1]*A[0][2];
c=A[2][0]*A[0][1]*A[1][2];
l=a+b+c;
a=A[2][0]*A[1][1]*A[0][2];
b=A[0][0]*A[2][1]*A[1][2];
c=A[1][0]*A[0][1]*A[2][2];
r=a+b+c;
res=l-r;
return res;
}


void cehk_circ_jesw(parm A,parm B,parm C,double *rad,parm *omega)
{int i;
double **MAT,a,b,A2,B2,C2,s,w;
parm S;
MAT=(double **)malloc(3*sizeof(double*));
for(i=0;i<3;i++)
	MAT[i]=(double *)malloc(3*sizeof(double));
MAT[0][0]=A.u;  MAT[0][1]=A.v;  MAT[0][2]=1.0;
MAT[1][0]=B.u;  MAT[1][1]=B.v;  MAT[1][2]=1.0;
MAT[2][0]=C.u;  MAT[2][1]=C.v;  MAT[2][2]=1.0;
a=wukb_dete_fudp(MAT);
A2=pojk_norm_jumt(A);
B2=pojk_norm_jumt(B);
C2=pojk_norm_jumt(C);
MAT[0][0]=A.u;  MAT[0][1]=A.v;  MAT[0][2]=A2;
MAT[1][0]=B.u;  MAT[1][1]=B.v;  MAT[1][2]=B2;
MAT[2][0]=C.u;  MAT[2][1]=C.v;  MAT[2][2]=C2;
b=wukb_dete_fudp(MAT);
MAT[0][0]=A2;	MAT[0][1]=A.v;  MAT[0][2]=1.0;
MAT[1][0]=B2;   MAT[1][1]=B.v;  MAT[1][2]=1.0;
MAT[2][0]=C2;   MAT[2][1]=C.v;  MAT[2][2]=1.0;
S.u=0.5*wukb_dete_fudp(MAT);
MAT[0][0]=A.u;  MAT[0][1]=A2;   MAT[0][2]=1.0;
MAT[1][0]=B.u;  MAT[1][1]=B2;   MAT[1][2]=1.0;
MAT[2][0]=C.u;  MAT[2][1]=C2;   MAT[2][2]=1.0;
S.v=0.5*wukb_dete_fudp(MAT);
for(i=0;i<3;i++)
	free(MAT[i]);
free(MAT);
omega->u=S.u/a;
omega->v=S.v/a;
s=sqrt(S.u*S.u+S.v*S.v);
w=(b/a)+((s*s)/(a*a));
*rad=sqrt(w);
}


int nefr_inte_sujc(point *A,double r,sphere *S1,sphere *S2)
{int i,suc=SUCCESS;
double s,d,temp,u,v,rad;
map_hyper M;
parm *a,omega;
point G;
vect3D nrm;
jags_find_mavk_nurp(A[0],A[1],A[2],&M);
a=(parm *)malloc(3*sizeof(parm));
for(i=0;i<3;i++)
	{suvd_prei_walj(M,A[i],&u,&v);
	a[i].u=u;
	a[i].v=v;
	}
cehk_circ_jesw(a[0],a[1],a[2],&rad,&omega);
free(a);
macn_imag_vuph(M,omega.u,omega.v,&G);
gotq_norm_bitg(A[0],A[1],A[2],&nrm);
s=wodt_dist_gilq(G,A[0]);
if(r<s)
	suc=FAILURE;
if(suc==SUCCESS)
	{temp=r*r-s*s;
	d=sqrt(temp);
	S1->zent.absi=G.absi+d*nrm.absi;
	S1->zent.ordo=G.ordo+d*nrm.ordo;
	S1->zent.cote=G.cote+d*nrm.cote;
	S1->rad=r;
	S2->zent.absi=G.absi-d*nrm.absi;
	S2->zent.ordo=G.ordo-d*nrm.ordo;
	S2->zent.cote=G.cote-d*nrm.cote;
	S2->rad=r;
	}
return suc;
}


void wusd_eval_jomk(c_arc3D C,double t,
double a,double b,point *P)
{double T;
T=t*b+(1.0-t)*a;
nehl_eval_segt(C,T,P);
}


void suzr_find_heqz_wikr(baryc2D q,c_arc3D C0,c_arc3D C1,
c_arc3D C2,point *sol)
{double a0,b0,a1,b1,a2,b2;
point term1,term2,term3,temp1,temp2,temp3;
qirp_inte_ligr(C0,&a0,&b0);
qirp_inte_ligr(C1,&a1,&b1);
qirp_inte_ligr(C2,&a2,&b2);
wusd_eval_jomk(C0,q.lambda2,a0,b0,&temp1);
wusd_eval_jomk(C2,1.0-q.lambda3,a2,b2,&temp2);
wusd_eval_jomk(C2,1.0,a2,b2,&temp3);
term1.absi=temp1.absi+temp2.absi-temp3.absi;
term1.ordo=temp1.ordo+temp2.ordo-temp3.ordo;
term1.cote=temp1.cote+temp2.cote-temp3.cote;
wusd_eval_jomk(C1,q.lambda3,a1,b1,&temp1);
wusd_eval_jomk(C0,1.0-q.lambda1,a0,b0,&temp2);
wusd_eval_jomk(C0,1.0,a0,b0,&temp3);
term2.absi=temp1.absi+temp2.absi-temp3.absi;
term2.ordo=temp1.ordo+temp2.ordo-temp3.ordo;
term2.cote=temp1.cote+temp2.cote-temp3.cote;
wusd_eval_jomk(C2,q.lambda1,a2,b2,&temp1);
wusd_eval_jomk(C1,1.0-q.lambda2,a1,b1,&temp2);
wusd_eval_jomk(C1,1.0,a1,b1,&temp3);
term3.absi=temp1.absi+temp2.absi-temp3.absi;
term3.ordo=temp1.ordo+temp2.ordo-temp3.ordo;
term3.cote=temp1.cote+temp2.cote-temp3.cote;
sol->absi=q.lambda1*term1.absi+q.lambda2*term2.absi+q.lambda3*term3.absi;
sol->ordo=q.lambda1*term1.ordo+q.lambda2*term2.ordo+q.lambda3*term3.ordo;
sol->cote=q.lambda1*term1.cote+q.lambda2*term2.cote+q.lambda3*term3.cote;
}


void kurl_find_vugz_sebt(baryc2D q,ns_curv C0,ns_curv C1,
ns_curv C2,parm *sol)
{parm term1,term2,term3,temp1,temp2,temp3;
kuqt_eval_webp(C0,q.lambda2,&temp1);
kuqt_eval_webp(C2,1.0-q.lambda3,&temp2);
kuqt_eval_webp(C2,1.0,&temp3);
term1.u=temp1.u+temp2.u-temp3.u;
term1.v=temp1.v+temp2.v-temp3.v;
kuqt_eval_webp(C1,q.lambda3,&temp1);
kuqt_eval_webp(C0,1.0-q.lambda1,&temp2);
kuqt_eval_webp(C0,1.0,&temp3);
term2.u=temp1.u+temp2.u-temp3.u;
term2.v=temp1.v+temp2.v-temp3.v;
kuqt_eval_webp(C2,q.lambda1,&temp1);
kuqt_eval_webp(C1,1.0-q.lambda2,&temp2);
kuqt_eval_webp(C1,1.0,&temp3);
term3.u=temp1.u+temp2.u-temp3.u;
term3.v=temp1.v+temp2.v-temp3.v;
sol->u=q.lambda1*term1.u+q.lambda2*term2.u+q.lambda3*term3.u;
sol->v=q.lambda1*term1.v+q.lambda2*term2.v+q.lambda3*term3.v;
}


void lekr_find_sujr_niqc(baryc2D q,c_curve cc,parm *sol)
{kurl_find_vugz_sebt(q,cc.nc[0],cc.nc[1],cc.nc[2],sol);
}

