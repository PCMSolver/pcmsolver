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
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "sas.h"
#include "partcas.h"



void rejm_coon_resd(c_curve cc,
double *T,double t,parm *X)
{double a,b,u,lambda;
parm temp;
a=T[0];
b=T[1];
lambda=t;
u=lambda*b+(1.0-lambda)*a;
mecg_eval_zupc(cc,u,&temp);
X->u=temp.u;
X->v=temp.v;
}



void rutj_coon_godh(c_curve cc,
double *T,double t,parm *X)
{double a,b,u,lambda;
parm temp;
a=T[1];
b=T[2];
lambda=t;
u=lambda*b+(1.0-lambda)*a;
mecg_eval_zupc(cc,u,&temp);
X->u=temp.u;
X->v=temp.v;
}



void fotp_coon_tojk(c_curve cc,
double *T,double t,parm *X)
{double a,b,u,lambda;
parm temp;
a=T[2];
b=T[3];
lambda=1.0-t;
u=lambda*b+(1.0-lambda)*a;
mecg_eval_zupc(cc,u,&temp);
X->u=temp.u;
X->v=temp.v;
}



void welr_coon_demp(c_curve cc,
double *T,double t,parm *X)
{double a,b,u,lambda;
parm temp;
a=T[3];
b=T[4];
lambda=1.0-t;
u=lambda*b+(1.0-lambda)*a;
mecg_eval_zupc(cc,u,&temp);
X->u=temp.u;
X->v=temp.v;
}


double fejl_find_gect_nomg(double t)
{double res;
res=t;
return res;
}


void fezd_coon_tekl(c_curve cc,
double u,double v,parm *sol)
{int nb_ct=4;
double F0u,F0v,F1u,F1v,A,B,C,*T,*PS;
parm alpha_u,gamma_u,delta_v,beta_v,alpha_0,alpha_1,gamma_0,gamma_1;
T =(double *)malloc((nb_ct+1)*sizeof(double));
PS=(double *)malloc((nb_ct+1)*sizeof(double));
jofv_dete_fatg(cc,T,PS);
rejm_coon_resd(cc,T,u,&alpha_u);
fotp_coon_tojk(cc,T,u,&gamma_u);
welr_coon_demp(cc,T,v,&delta_v);
rutj_coon_godh(cc,T,v,&beta_v);
rejm_coon_resd(cc,T,0.0,&alpha_0);
rejm_coon_resd(cc,T,1.0,&alpha_1);
fotp_coon_tojk(cc,T,0.0,&gamma_0);
fotp_coon_tojk(cc,T,1.0,&gamma_1);
free(T);
free(PS);
F1u=fejl_find_gect_nomg(u);
F0u=1.0-F1u;
F1v=fejl_find_gect_nomg(v);
F0v=1.0-F1v;
A=alpha_u.u*F0v+gamma_u.u*F1v;
B=-delta_v.u+alpha_0.u*F0v+gamma_0.u*F1v;
C=-beta_v.u+alpha_1.u*F0v+gamma_1.u*F1v;
sol->u=A-F0u*B-F1u*C;
A=alpha_u.v*F0v+gamma_u.v*F1v;
B=-delta_v.v+alpha_0.v*F0v+gamma_0.v*F1v;
C=-beta_v.v+alpha_1.v*F0v+gamma_1.v*F1v;
sol->v=A-F0u*B-F1u*C;
}


void wudj_allo_worv(int nnd,int nel,int ned,quadrangulation *quad)
{quad->knot=(parm *)malloc(nnd*sizeof(parm));
quad->elem=(efajor *)malloc(nel*sizeof(efajor));
quad->kt=(kt *)malloc(ned*sizeof(kt));
quad->zt=(double *)malloc(nnd*sizeof(double));
quad->flag=(int *)malloc(nnd*sizeof(int));
quad->type=(int *)malloc(nnd*sizeof(int));
}


void wacm_dest_wokn(quadrangulation *quad)
{free(quad->knot);
free(quad->elem);
free(quad->kt);
free(quad->zt);
free(quad->flag);
free(quad->type);
}


void gedp_inte_vewq(manif_ro unt,quadrangulation quad,
trmsrf surf,c_curve c_loc,
double rd,double gr,double bl,
manif_tl *msh,rgb_lk *col)
{int NND,NEL,nnd,nel,i,*map;
int n1,n2,n3;
double u,v;
parm pr;
point img;
NND=msh->n_grs;
NEL=msh->e_grs;
nnd=unt.n_grs;
nel=unt.e_grs;
map=(int *)malloc(nnd*sizeof(int));
for(i=0;i<nnd;i++)
	{u=unt.knot[i].u;
	v=unt.knot[i].v;
	fezd_coon_tekl(c_loc,u,v,&pr);
	wolf_eval_murg(surf,pr.u,pr.v,&img);
	getf_find_rogc_todj(img,&msh->knot[NND]);
	map[i]=NND;
	NND++;
	}
for(i=0;i<nel;i++)
	{n1=unt.entity[i].frvrt;
	n2=unt.entity[i].scvrt;
	n3=unt.entity[i].thvrt;
	msh->entity[NEL].frvrt=map[n1];
	msh->entity[NEL].scvrt=map[n2];
	msh->entity[NEL].thvrt=map[n3];
	col[NEL].red=rd;
	col[NEL].green=gr;
	col[NEL].blue=bl;
	NEL++;
	}
free(map);
msh->n_grs=NND;
msh->e_grs=NEL;
}


void sotp_atom_kegt(atom *A,int nb)
{int N=5,nnd_loc,nel_loc,k_el;
int p,q,i,s,nel,hemi[2];
double rd,gr,bl;
prop_ccurve pcc;
trmsrf surf;
quadrangulation quad;
c_curve c_loc;
rgb_lk *col;
manif_ro unt;
manif_tl msh;
nnd_loc    =N*N;
nel_loc    =2*(N-1)*(N-1);
unt.knot   =(parm *)malloc(nnd_loc*sizeof(parm));
unt.entity=(telolf *)malloc(nel_loc*sizeof(telolf));
lenq_simp_socg(N,N,&unt);
hemi[0]=NORTH_HEMI;
hemi[1]=SOUTH_HEMI;
pcc.N=4;
for(i=0;i<4;i++)
	{pcc.k[i]=5;
	pcc.n[i]=10;
	}
pcc.nle=4;
pcc.nca=4;
pcc.nnc=4;
homd_allo_tevf(pcc,&surf.cc);
wudj_allo_worv(8,5,20,&quad);
homd_allo_tevf(pcc,&c_loc);
msh.n_grs=0;
msh.e_grs=0;
msh.knot   =(point *)malloc(10*nnd_loc*nb*sizeof(point));
msh.entity=(telolf *)malloc(10*nel_loc*nb*sizeof(telolf));
col        =(rgb_lk *)malloc(10*nel_loc*nb*sizeof(rgb_lk));
k_el=0;
for(q=0;q<nb;q++)
for(p=0;p<2;p++)
	{kucw_atom_nocd(hemi[p],A[q],&surf);
	rahv_atom_zelf(surf,&quad);
	corm_fill_ruqt(&quad);
	nel=quad.e_grs;
	for(s=0;s<nel;s++)
		{gojw_quad_wuln(hemi[p],surf,quad,s,&c_loc);
		tesr_colo_donr(k_el,&rd,&gr,&bl);
		k_el++;
		gedp_inte_vewq(unt,quad,surf,c_loc,rd,gr,bl,&msh,col);
		}
	}
wacm_dest_wokn(&quad);
wosn_dest_jomw(pcc,&surf.cc);
wosn_dest_jomw(pcc,&c_loc);
free(unt.entity);
free(unt.knot);
free(msh.entity);
free(msh.knot);
free(col);
}


int varm_dete_qumz(atom *S,int N,double probe,int *list)
{int i,j,nb,ind,ts;
double D,R_i,R_j;
nb=0;
for(i=0;i<N;i++)
	{ind=1;
	for(j=0;j<N;j++)if(i!=j)
		{R_i= S[i].rad+probe;
		R_j = S[j].rad+probe;
		D   = R_i+R_j;
		ts=gect_tole_husn(S[i].zent,S[j].zent,D);
		if(ts==1)
			{ind=2;
			break;
			}
		}
	if(ind==1)
		{list[nb]=i;
		nb++;
		}
	}
return nb;
}

