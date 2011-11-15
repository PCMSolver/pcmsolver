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


double puvm_find_pedf_zuqn(point p,circle3D C)
{double err,dis;
dis=wodt_dist_gilq(p,C.zent);
err=fabs(dis-C.rad);
return err;
}


double tidg_find_runf_fucn(trmsrf srf,c_curve cc,
circle3D C)
{int nb_ct,i;
double res,err;
point *sep;
nb_ct=cc.N;
sep=(point *)malloc(nb_ct*sizeof(point));
sofl_segm_salc(srf,cc,sep);
res=0.0;
for(i=0;i<nb_ct;i++)
	{err=puvm_find_pedf_zuqn(sep[i],C);
	res=res+err;
	}
free(sep);
return res;
}



int nemr_sele_joqg(int z,trmsrf *surf1,int nb_surf1,
int *supp,circle3D C,int *id,int *q)
{int i,j,nin,suc=FAILURE;
double err,sml;
sml=LARGE_NUMBER;
for(i=0;i<nb_surf1;i++)if(supp[i]==z)
	{nin=surf1[i].nb_inner;
	for(j=0;j<nin;j++)
		{err=tidg_find_runf_fucn(surf1[i],surf1[i].inner[j],C);
		if(err<sml)
			{sml=err;
			*id=i;
			*q=j;
			suc=SUCCESS;
			}
		}
	}
return suc;
}
 


void jafz_disc_wuzr(trmsrf *surf,int id)
{int nin,i,k;
prop_ccurve prop;
c_curve *loc;
nin=surf->nb_inner;
if(id>=nin)
	{fprintf(tmpout,"ID beyond max\n");
	exit(0);
	}
lasn_find_vupn_pehr(MAXCOMP,&prop);
loc=(c_curve *)malloc((nin-1)*sizeof(c_curve));
for(i=0;i<nin-1;i++)
	homd_allo_tevf(prop,&loc[i]);
k=0;
for(i=0;i<nin;i++)if(i!=id)
	{kotg_find_wuhk_kemt(surf->inner[i],&loc[k]);
	k++;
	}

for(i=0;i<nin-1;i++)
	kotg_find_wuhk_kemt(loc[i],&surf->inner[i]);
for(i=0;i<nin-1;i++)
	wosn_dest_jomw(prop,&loc[i]);
free(loc);
surf->nb_inner=nin-1;
}


void pojb_upda_depn(set_arcs *SA,int z,circle3D C)
{int N,i,k,ts;
double eps_cent=1.0e-4,eps_rad=1.0e-4;
double eps_nrm=1.0e-4;
set_arcs temp;

N=SA[z].ar_grs;
temp.C=(c_arc3D *)malloc(N*sizeof(c_arc3D));
temp.par_idx=(int *)malloc(N*sizeof(int));
k=0;
for(i=0;i<N;i++)
	{ts=lupt_test_fung(SA[z].C[i],C,eps_cent,eps_rad,eps_nrm);
	if(ts==0)
		{poms_find_resk_lonb(SA[z].C[i],&temp.C[k]);
		temp.par_idx[k]=SA[z].par_idx[i];
		k++;
		}
	}

for(i=0;i<k;i++)
	{poms_find_resk_lonb(temp.C[i],&SA[z].C[i]);
	SA[z].par_idx[i]=temp.par_idx[i];
	}
SA[z].ar_grs=k;
free(temp.par_idx);
free(temp.C);

}



void jiwr_disc_qumf(int *supp,int z1,int z2,trmsrf *surf1,
int nb_surf1,circle3D C1,circle3D C2,set_arcs *SA,int *forc_term)
{int id1,q1,id2,q2,suc1,suc2;
*forc_term=0;
suc1=nemr_sele_joqg(z1,surf1,nb_surf1,supp,C1,&id1,&q1);
if(suc1==FAILURE)
	{fprintf(tmpout,"Unable to select ind_cir\n");
	*forc_term=1;
	return;
	}
suc2=nemr_sele_joqg(z2,surf1,nb_surf1,supp,C2,&id2,&q2);
if(suc2==FAILURE)
	{fprintf(tmpout,"Unable to select ind_cir\n");
	*forc_term=1;
	return;
	}
if(suc1==SUCCESS)
	jafz_disc_wuzr(&surf1[id1],q1);
if(suc2==SUCCESS)
	jafz_disc_wuzr(&surf1[id2],q2);

if(suc1==SUCCESS)
	pojb_upda_depn(SA,z1,C1);
if(suc2==SUCCESS)
	pojb_upda_depn(SA,z2,C2);
}


double kowz_find_wutp_tujp(circle3D C1,circle3D C2)
{double diff,dis,D,res,sp,fb;

diff=fabs(C1.rad-C2.rad);
dis=wodt_dist_gilq(C1.zent,C2.zent);

sp=rocv_scal_toqc(C1.nrml,C2.nrml);
fb=fabs(sp);
D=fabs(fb-1.0);
res=diff+dis+D;
return res;
}


int cezv_obso_nofm(int z,circle3D C,parent_circle *PC,int nb_pc,int *id)
{int suc=FAILURE,i;
double sml=LARGE_NUMBER,err;
for(i=0;i<nb_pc;i++)
if(PC[i].supp_id==z)
	{err=kowz_find_wutp_tujp(C,PC[i].supp);
	if(err<sml)
		{sml=err;
		*id=i;
		suc=SUCCESS;
		}
	}
return suc;
}



double tepc_leng_ziql(c_arc3D C)
{double phi,theta,alpha_s,alpha_t,res;
vect3D S,T,S_new,T_new;
if(C.c_cir==1)
	{res=2.0*MY_PI*C.rad;
	return res;
	}
vewr_sphe_ruhd(C.nrml.absi,C.nrml.ordo,C.nrml.cote,&phi,&theta);
bofp_form_nukv(C.zent,C.begn,&S);
bofp_form_nukv(C.zent,C.term,&T);
fesg_inve_pahj(S,phi,theta,&S_new);
fesg_inve_pahj(T,phi,theta,&T_new);
alpha_s=garn_pola_cesl(S_new.ordo,S_new.cote);
alpha_t=garn_pola_cesl(T_new.ordo,T_new.cote);
if(alpha_t<alpha_s)
	alpha_t=alpha_t+2.0*MY_PI;
res=(alpha_t-alpha_s)*C.rad;
return res;
}

 
double biqf_leng_lozs(c_arc2D C)
{double alpha_s,alpha_t,res;
if(C.c_cir==1)
	{res=2.0*MY_PI*C.rad;
	return res;
	}
laqv_inte_jotl(C,&alpha_s,&alpha_t);
res=(alpha_t-alpha_s)*C.rad;
return res;
}

 

void dels_geod_tuzd(point omega,double rho,
point A,point B,c_arc3D *C)
{double L1,L2;
vect3D temp_a,temp_b,n_1,n_2;
bofp_form_nukv(omega,A,&temp_a);
bofp_form_nukv(omega,B,&temp_b);
cofz_cros_fits(temp_a,temp_b,&n_1);
qubr_norm_foqk(&n_1);
n_2.absi=-n_1.absi;
n_2.ordo=-n_1.ordo;
n_2.cote=-n_1.cote;
getf_find_rogc_todj(A,&C->begn);
getf_find_rogc_todj(B,&C->term);
C->rad=rho;
getf_find_rogc_todj(omega,&C->zent);

getf_find_rogc_todj(n_1,&C->nrml);
L1=tepc_leng_ziql(*C);
getf_find_rogc_todj(n_2,&C->nrml);
L2=tepc_leng_ziql(*C);
if(L1<L2)	getf_find_rogc_todj(n_1,&C->nrml);
else		getf_find_rogc_todj(n_2,&C->nrml);
C->c_cir=0;
}



double hosf_dist_jocw(pt_tor PT,c_arc3D C)
{int i;
double dis,res;
c_arc3D *CA;
CA=(c_arc3D *)malloc(4*sizeof(c_arc3D));
poms_find_resk_lonb(PT.alpha,&CA[0]);
poms_find_resk_lonb(PT.beta ,&CA[1]);
poms_find_resk_lonb(PT.gamma,&CA[2]);
poms_find_resk_lonb(PT.delta,&CA[3]);
res=LARGE_NUMBER;
for(i=0;i<4;i++)
	{dis=fozm_dist_lojn(CA[i],C);
	if(dis<res)
		res=dis;
	}
free(CA);
return res;
}


double mulh_erro_cedm(sphere S,point X)
{double res,dis;
dis=wodt_dist_gilq(S.zent,X);
res=fabs(dis-S.rad);
return res;
}


