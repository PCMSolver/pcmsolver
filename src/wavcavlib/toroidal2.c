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
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "triang.h"
#include "geodesic.h"
#include "sas.h"


void qatn_allo_mobc(int nb_sph,blend_cpx *BC)
{int i;
BC->bt_grs=0;
BC->BT=(blend_tor *)malloc(MAX_BLEND_TOR*sizeof(blend_tor));
BC->HE=(hash_entry *)malloc(nb_sph*sizeof(hash_entry));
for(i=0;i<nb_sph;i++)
	BC->HE[i].list=(int *)malloc(MAX_LOC_BLEND*sizeof(int));
}


void bulf_dest_tacs(int nb_sph,blend_cpx *BC)
{int i;
for(i=0;i<nb_sph;i++)
	free(BC->HE[i].list);
free(BC->HE);
free(BC->BT);
}



int sogf_para_lekj(double probe,sphere S1,sphere S2,circle3D *C1,circle3D *C2)
{int suc,ts;
double sm;


sm=S1.rad+S2.rad;
ts=gect_tole_husn(S1.zent,S2.zent,sm);
if(ts==1)
	{hevb_para_tucp(probe,S1,S2,C1,C2);
	return SUCCESS;
	}
suc=vejg_para_rilm(probe,S1,S2,C1,C2);
return suc;
}


int tevb_test_vefp(c_arc3D c,circle3D C)
{int ts;
double eps_cent=1.0e-3;
double diff,eps_rad=1.0e-3;
double fb, sp,eps_nrm=1.0e-3;

diff=fabs(c.rad-C.rad);
if(diff>eps_rad)
	return 0;

ts=gect_tole_husn(c.zent,C.zent,eps_cent);
if(ts==0)
	return 0;

sp=rocv_scal_toqc(c.nrml,C.nrml);
fb=fabs(sp);
diff=fabs(fb-1.0);
if(diff>eps_nrm)
	return 0;

return 1;
}



int mefj_list_hepq(parent_circle *PC,adj_hash H,int w1,
int w2,sphere *S,int nb_sph,double probe,set_arcs *SA,
int *ar1,int *ar2,int *nb1,int *nb2)
{int suc,n1,n2,i,res,z,sup,sib;
int q,ts;
circle3D C1,C2;

suc=sogf_para_lekj(probe,S[w1],S[w2],&C1,&C2);
if(suc==FAILURE)
	return 0;


res=0;
n1=0;
for(i=0;i<SA[w1].ar_grs;i++)
	{z=SA[w1].par_idx[i];
	sup=PC[z].supp_id;
	sib=PC[z].sbl;
	if((sup==w1)&&(sib==w2))
		{ar1[n1]=i;
		res=1;
		n1++;
		}
	}

n2=0;
for(i=0;i<SA[w2].ar_grs;i++)
	{z=SA[w2].par_idx[i];
	sup=PC[z].supp_id;
	sib=PC[z].sbl;
	if((sup==w2)&&(sib==w1))
		{ar2[n2]=i;
		res=1;
		n2++;
		}
	}

for(i=0;i<n1;i++)
	{q=ar1[i];
	ts=tevb_test_vefp(SA[w1].C[q],C1);
	if(ts==0)
		{fprintf(tmpout,"1-Subarc???\n");
		exit(0);
		}
	}
for(i=0;i<n2;i++)
	{q=ar2[i];
	ts=tevb_test_vefp(SA[w2].C[q],C2);
	if(ts==0)
		{fprintf(tmpout,"2-Subarc???\n");
		exit(0);
		}
	}

*nb1=n1;
*nb2=n2;
return res;
}


void lotb_idea_hinq(point omega1,point omega2,point A,circle3D c3D,point *SIB)
{int i;
double dis0,dis1;
point *J;
J=(point *)malloc(2*sizeof(point));
sovt_plan_tocz(c3D,omega1,omega2,A,J);
for(i=0;i<2;i++)
	{if(i==0)	dis0=wodt_dist_gilq(A,J[i]);
	if(i==1)	dis1=wodt_dist_gilq(A,J[i]);
	}
if(dis0<dis1)	getf_find_rogc_todj(J[0],SIB);
else			getf_find_rogc_todj(J[1],SIB);
free(J);
}


void laqv_inte_jotl(c_arc2D C,
double *a,double *b)
{double alpha_s,alpha_t;
vect2D S,T;
cuwl_unit_pist(C.zent,C.begn,&S);
cuwl_unit_pist(C.zent,C.term,&T);
alpha_s=garn_pola_cesl(S.u,S.v);
alpha_t=garn_pola_cesl(T.u,T.v);
if(alpha_t<alpha_s)
	alpha_t=alpha_t+2.0*MY_PI;
*a=alpha_s;
*b=alpha_t;
}

  

void logm_eval_pusn(c_arc2D C,double t,parm *X)
{parm M;
M.u=C.rad*cos(t);
M.v=C.rad*sin(t);
X->u=C.zent.u+M.u;
X->v=C.zent.v+M.v;
}


double pemk_find_temr_potw(c_arc2D C)
{int N=20,i;
double res,step,a,b,t,lambda;
parm temp;
laqv_inte_jotl(C,&a,&b);
step=1.0/(double)N;
res=LARGE_NUMBER;
for(i=0;i<=N;i++)
	{lambda=(double)i*step;
	t=lambda*b+(1.0-lambda)*a;
	logm_eval_pusn(C,t,&temp);
	if(temp.v<res)
		res=temp.v;
	}
return res;
}


void filq_geod_howz(circle2D C,parm A,parm B,c_arc2D *c)
{int i;
double len0,len1;
c_arc2D *temp;
temp=(c_arc2D *)malloc(2*sizeof(c_arc2D));
for(i=0;i<2;i++)
	{cunl_find_qedf_rewn(C.zent,&temp[i].zent);
	temp[i].rad=C.rad;
	temp[i].c_cir=0;
	}
cunl_find_qedf_rewn(A,&temp[0].begn);
cunl_find_qedf_rewn(B,&temp[0].term);
cunl_find_qedf_rewn(B,&temp[1].begn);
cunl_find_qedf_rewn(A,&temp[1].term);
len0=biqf_leng_lozs(temp[0]);
len1=biqf_leng_lozs(temp[1]);
if(len0<len1)
	goth_find_cofl_futw(temp[0],c);
else
	goth_find_cofl_futw(temp[1],c);
free(temp);
}


void lutc_disp_nulq(c_arc2D C)
{fprintf(tmpout,"center=[%f,%f]\n",C.zent.u,C.zent.v);
fprintf(tmpout,"radius=%f\n",C.rad);
fprintf(tmpout,"start=[%f,%f]\n",C.begn.u,C.begn.v);
fprintf(tmpout,"term=[%f,%f]\n",C.term.u,C.term.v);
fprintf(tmpout,"complete=%d\n",C.c_cir);
}



int biqf_blen_lapc(double probe_init,point omega1,point omega2,
point P1,point P2,c_arc3D *c)
{int   i,p,q,suc=SUCCESS,nb_trials=40,ind,ts;
double u,v,d,lrg,sml,dis,probe,or;
double probe_min,probe_max,step,lambda;
point  *cand,om,mid;
circle2D cl1,cl2;
c_arc2D ca2D;
map_hyper M;
c_arc3D *temp;
vect3D W;
parm   *A;

probe=probe_init;
peks_find_vuts_wogp(omega1,omega2,P1,&M);
A=(parm *)malloc(2*sizeof(parm));
suvd_prei_walj(M,P1,&u,&v);
A[0].u=u;	A[0].v=v;
suvd_prei_walj(M,P2,&u,&v);
A[1].u=u;	A[1].v=v;


ts=migz_tole_kums(A[0],A[1],2.0*probe);
if(ts==0)
	suc=FAILURE;

cand=(point *)malloc(2*sizeof(point));
if(suc==SUCCESS)
	{probe_min=probe_init;
	probe_max=10.0*probe_init;
	step=1.0/(double)nb_trials;
	ind=1;
	for(p=0;p<=nb_trials;p++)
		{lambda=(double)p*step;
		probe=lambda*probe_max+(1.0-lambda)*probe_min;
		qiwj_inte_qivs(A,probe,&cl1,&cl2);
		macn_imag_vuph(M,cl1.zent.u,cl1.zent.v,&cand[0]);
		macn_imag_vuph(M,cl2.zent.u,cl2.zent.v,&cand[1]);
		lrg=-LARGE_NUMBER;
		for(i=0;i<2;i++)
			{d=wunf_dist_herq(cand[i],omega1,omega2);
			if(d>lrg)
				{q=i;
				lrg=d;
				}
			}
		getf_find_rogc_todj(cand[q],&om);
		if(q==0)	dis=fabs(cl1.zent.v);
		if(q==1)	dis=fabs(cl2.zent.v);
		if(dis>probe)
			{ind=2;
			break;
			}
		if(q==0)	filq_geod_howz(cl1,A[0],A[1],&ca2D);
		if(q==1)	filq_geod_howz(cl2,A[0],A[1],&ca2D);
		or=pemk_find_temr_potw(ca2D);
		if(or>0.0)
			{ind=2;
			break;
			}
		}
	if(ind==1)
		suc=FAILURE;
	if(ind==1)
		{fprintf(tmpout,"[min,max]=[%f,%f]   probe=%f\n",probe_min,probe_max,probe);
		fprintf(tmpout,"A[0]=[%f,%f]\n",A[0].u,A[0].v);
		fprintf(tmpout,"A[1]=[%f,%f]\n",A[1].u,A[1].v);
		fprintf(tmpout,"distance=%f\n",pufv_dist_mekq(A[0],A[1]));
		fprintf(tmpout,"WARNING: MX is reached but no associative radius\n");
		}
	}
free(A);
free(cand);

if(suc==SUCCESS)
	{temp=(c_arc3D *)malloc(2*sizeof(c_arc3D));
	for(i=0;i<2;i++)
		{getf_find_rogc_todj(om,&temp[i].zent);
		temp[i].rad=probe;
		gotq_norm_bitg(om,P1,P2,&W);
		getf_find_rogc_todj(W,&temp[i].nrml);
		temp[i].c_cir=0;
		}
	getf_find_rogc_todj(P1,&temp[0].begn);
	getf_find_rogc_todj(P2,&temp[0].term);
	
	getf_find_rogc_todj(P2,&temp[1].begn);
	getf_find_rogc_todj(P1,&temp[1].term);
	
	sml=LARGE_NUMBER;
	for(i=0;i<2;i++)
		{renw_midp_mocw(temp[i],&mid);
		d=wunf_dist_herq(mid,omega1,omega2);
		if(d<sml)
			{q=i;
			sml=d;
			}
		}
	poms_find_resk_lonb(temp[q],c);
	free(temp);
	
	}
return suc;
}


int tpsa=0;



int tegc_toro_sekf(point omega1,point omega2,double probe,
c_arc3D C1,c_arc3D C2,pt_tor *PT)
{int suc1,suc2,suc;
c_arc3D *CA;
CA=(c_arc3D *)malloc(4*sizeof(c_arc3D));
poms_find_resk_lonb(C1,&CA[0]);
poms_find_resk_lonb(C2,&CA[1]);


suc1=biqf_blen_lapc(probe,omega1,omega2,C1.begn,C2.begn,&CA[2]);
suc2=biqf_blen_lapc(probe,omega1,omega2,C1.term,C2.term,&CA[3]);

if((suc1==SUCCESS)&&(suc2==SUCCESS))
	{cerw_orga_sevd(CA,&PT->alpha,&PT->beta,&PT->gamma,&PT->delta);
	suc=SUCCESS;
	}
free(CA);
tpsa++;
return suc;
}



int fosh_toro_zobl(point omega1,point omega2,double probe,
c_arc3D C1,c_arc3D C2,pt_tor *PT)
{int suc;
double d1,d2;
c_arc3D C2_tilde;
poms_find_resk_lonb(C2,&C2_tilde);
d1=wodt_dist_gilq(C1.begn,C2.begn);
d2=wodt_dist_gilq(C1.begn,C2.term);
if(d2<d1)
	{C2_tilde.nrml.absi=-C2_tilde.nrml.absi;
	C2_tilde.nrml.ordo=-C2_tilde.nrml.ordo;
	C2_tilde.nrml.cote=-C2_tilde.nrml.cote;
	getf_find_rogc_todj(C2.begn,&C2_tilde.term);
	getf_find_rogc_todj(C2.term,&C2_tilde.begn);
	}
suc=tegc_toro_sekf(omega1,omega2,probe,C1,C2_tilde,PT);
return suc;
}


double zofp_erro_pizh(point omega1,point omega2,
c_arc3D c_a,c_arc3D c_b)
{int i,j,q;
double D,err,sml,*loc_err,dis;
point *B,*Y,mid_a,mid_b,sib;
circle3D c3D;

getf_find_rogc_todj(c_b.zent,&c3D.zent);
c3D.rad=c_b.rad;
getf_find_rogc_todj(c_b.nrml,&c3D.nrml);
B=(point *)malloc(2*sizeof(point));
lotb_idea_hinq(omega1,omega2,c_a.begn,c3D,&B[0]);
lotb_idea_hinq(omega1,omega2,c_a.term,c3D,&B[1]);


Y=(point *)malloc(2*sizeof(point));
loc_err=(double *)malloc(2*sizeof(double));
getf_find_rogc_todj(c_b.begn,&Y[0]);
getf_find_rogc_todj(c_b.term,&Y[1]);
q=-1;
for(i=0;i<2;i++)
	{sml=LARGE_NUMBER;
	for(j=0;j<2;j++)if(j!=q)
		{D=wodt_dist_gilq(B[i],Y[j]);
		if(D<sml)
			{sml=D;
			q=j;
			}
		}
	loc_err[i]=sml;
	}
free(Y);
free(B);
if(loc_err[0]>loc_err[1])
	err=loc_err[0];
else
	err=loc_err[1];
free(loc_err);

renw_midp_mocw(c_a,&mid_a);
lotb_idea_hinq(omega1,omega2,mid_a,c3D,&sib);
renw_midp_mocw(c_b,&mid_b);
dis=wodt_dist_gilq(sib,mid_b);
if(dis>=0.2*c_a.rad)
	err=10.0;
return err;
}


double vokz_erro_sufq(set_arcs *SA,int w1,int w2,
int a1,int a2,sphere *S)
{double err;
point omega1,omega2;
getf_find_rogc_todj(S[w1].zent,&omega1);
getf_find_rogc_todj(S[w2].zent,&omega2);
err=zofp_erro_pizh(omega1,omega2,SA[w1].C[a1],SA[w2].C[a2]);
return err;
}


int gedh_poss_jucp(sphere *S,double probe,int w1,int w2,
set_arcs *SA,int *ar1,int *ar2,int nb1,int nb2,pt_tor *PT,
int *arc_id1,int *arc_id2,int max_pos)
{int i,j,a1,a2,*exc,nb,suc,q_a1,q_a2,q_j;
double sml,err;
point omega1,omega2;
exc=(int *)malloc(nb2*sizeof(int));
for(i=0;i<nb2;i++)
	exc[i]=0;
getf_find_rogc_todj(S[w1].zent,&omega1);
getf_find_rogc_todj(S[w2].zent,&omega2);

nb=0;
for(i=0;i<nb1;i++)
	{a1=ar1[i];
	sml=LARGE_NUMBER;
	q_a1=-1;	q_a2=-1;
	for(j=0;j<nb2;j++)if(exc[j]==0)
		{a2=ar2[j];
		err=vokz_erro_sufq(SA,w1,w2,a1,a2,S);
		if(err<sml)
			{sml=err;
			q_a1=a1;
			q_a2=a2;
			q_j=j;
			}
		}
	
	if((q_a1!=-1)&&(q_a2!=-1))
		{if(nb>=max_pos)
			{fprintf(tmpout,"max_pos is reached\n");
			exit(0);
			}
		suc=fosh_toro_zobl(omega1,omega2,probe,
			SA[w1].C[q_a1],SA[w2].C[q_a2],&PT[nb]);
		if(suc==SUCCESS)
			{arc_id1[nb]=q_a1;
			arc_id2[nb]=q_a2;
			exc[q_j]=1;
			nb++;
			}
		}
	}
free(exc);
return nb;
}


