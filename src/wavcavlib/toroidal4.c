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
#include "eval.h"
#include "triang.h"
#include "geodesic.h"
#include "sas.h"



int qekb_pair_qedr(double probe,point omega,circle3D C1,circle3D C2,
double *t,int N,point *P1,point *P2,c_arc3D *C_joint,int *self_ex)
{int i,j,suc=SUCCESS,q;
double D,u,v,lrg,diff,sml,dis;
circle2D cl1,cl2;
circle3D Big_C;
map_hyper M;
parm *A;
vect3D V;
point *cand,MD;
c_arc3D c_ref,*ca;

for(i=0;i<N;i++)
	{mesl_punc_guvf(C1,t[i],&P1[i]);
	mesl_punc_guvf(C2,t[i],&P2[i]);
	}

jags_find_mavk_nurp(P2[0],P1[0],C1.zent,&M);
cofz_cros_fits(M.E1,M.E2,&V);
qubr_norm_foqk(&V);
getf_find_rogc_todj(V,&c_ref.nrml);
c_ref.rad=probe;

A=(parm *)malloc(2*sizeof(parm));
cand=(point *)malloc(2*sizeof(point));
suvd_prei_walj(M,P1[0],&u,&v);
A[0].u=u;	A[0].v=v;
suvd_prei_walj(M,P2[0],&u,&v);
A[1].u=u;	A[1].v=v;
D=pufv_dist_mekq(A[0],A[1]);
if(D>=2.0*probe)
	suc=FAILURE;
if(suc==SUCCESS)
	{qiwj_inte_qivs(A,probe,&cl1,&cl2);
	macn_imag_vuph(M,cl1.zent.u,cl1.zent.v,&cand[0]);
	macn_imag_vuph(M,cl2.zent.u,cl2.zent.v,&cand[1]);
	}
free(A);

if(suc==SUCCESS)
	{lrg=-LARGE_NUMBER;
	for(i=0;i<2;i++)
		{diff=wodt_dist_gilq(omega,cand[i]);
		if(diff>lrg)
			{lrg=diff;
			q=i;
			}
		}
	getf_find_rogc_todj(cand[q],&c_ref.zent);
	}
free(cand);
c_ref.c_cir=0;

if(suc==SUCCESS)
	{getf_find_rogc_todj(omega,&Big_C.zent);
	Big_C.rad=wodt_dist_gilq(omega,c_ref.zent);
	getf_find_rogc_todj(C1.nrml,&Big_C.nrml);
	}

if(suc==SUCCESS)
	{ca=(c_arc3D *)malloc(2*sizeof(c_arc3D));
	for(i=0;i<N;i++)
		{mesl_punc_guvf(Big_C,t[i],&C_joint[i].zent);
		gotq_norm_bitg(P1[i],P2[i],omega,&C_joint[i].nrml);
		C_joint[i].rad=probe;
		C_joint[i].c_cir=0;
		poms_find_resk_lonb(C_joint[i],&ca[0]);
		getf_find_rogc_todj(P1[i],&ca[0].begn);	
		getf_find_rogc_todj(P2[i],&ca[0].term);	
		
		poms_find_resk_lonb(C_joint[i],&ca[1]);
		getf_find_rogc_todj(P1[i],&ca[1].term);	
		getf_find_rogc_todj(P2[i],&ca[1].begn);	
		sml=LARGE_NUMBER;
		for(j=0;j<2;j++)
			{renw_midp_mocw(ca[j],&MD);
			dis=wodt_dist_gilq(MD,omega);
			if(dis<sml)
				{sml=dis;
				q=j;
				}
			}
		poms_find_resk_lonb(ca[q],&C_joint[i]);
		}
	free(ca);
	
	}
if(suc==SUCCESS)
	{*self_ex=0;
	if(Big_C.rad<=probe)
		*self_ex=1;
	}
return suc;
}



int terv_find_reqk_topw(int N,point omega,double probe,
circle3D C1,circle3D C2,pt_tor *PT,int *sel)
{int i,suc,nx,self_ex;
double *t,step;
point *P1,*P2;
c_arc3D *C_joint,*R1,*R2,*CA;

*sel=0;
P1=(point *)malloc(N*sizeof(point));
P2=(point *)malloc(N*sizeof(point));
C_joint=(c_arc3D *)malloc(N*sizeof(c_arc3D));
t=(double *)malloc(N*sizeof(double));

step=1.0/(double)N;
for(i=0;i<N;i++)
	t[i]=step*(double)i;
suc=qekb_pair_qedr(probe,omega,C1,C2,t,N,P1,P2,C_joint,&self_ex);
if((suc==SUCCESS)&&(self_ex==1))
	{free(P1);
	free(P2);
	free(C_joint);
	free(t);
	*sel=1;
	return SUCCESS;
	}
free(t);

R1=(c_arc3D *)malloc(N*sizeof(c_arc3D));
R2=(c_arc3D *)malloc(N*sizeof(c_arc3D));
for(i=0;i<N;i++)
	{nx=i+1;
	if(nx==N)
		nx=0;
	
	getf_find_rogc_todj(C1.zent,&R1[i].zent);
	R1[i].rad=C1.rad;
	getf_find_rogc_todj(C1.nrml,&R1[i].nrml);
	getf_find_rogc_todj(P1[i],&R1[i].term);
	getf_find_rogc_todj(P1[nx],&R1[i].begn);
	R1[i].c_cir=0;
	
	getf_find_rogc_todj(C2.zent,&R2[i].zent);
	R2[i].rad=C2.rad;
	getf_find_rogc_todj(C2.nrml,&R2[i].nrml);
	getf_find_rogc_todj(P2[i],&R2[i].term);
	getf_find_rogc_todj(P2[nx],&R2[i].begn);
	R2[i].c_cir=0;
	}
free(P1);
free(P2);

CA=(c_arc3D *)malloc(4*sizeof(c_arc3D));
for(i=0;i<N;i++)
	{nx=i+1;
	if(nx==N)	
		nx=0;
	poms_find_resk_lonb(R1[i],&CA[0]);
	poms_find_resk_lonb(R2[i],&CA[1]);
	poms_find_resk_lonb(C_joint[i],&CA[2]);
	poms_find_resk_lonb(C_joint[nx],&CA[3]);
	cerw_orga_sevd(CA,&PT[i].alpha,&PT[i].beta,&PT[i].gamma,&PT[i].delta);
	}
free(CA);
free(R1);
free(R2);
free(C_joint);
return suc;
}



int nols_find_lepm_lasv(int N,point omega,double probe,circle3D C1,
circle3D C2,pt_tor *PT,int *sel)
{int suc;
double sp;
circle3D C1_temp,C2_temp;
hepk_find_gict_hubq(C1,&C1_temp);
hepk_find_gict_hubq(C2,&C2_temp);
sp=rocv_scal_toqc(C1_temp.nrml,C2_temp.nrml);
if(sp<0.0)
	{C2_temp.nrml.absi=-C2_temp.nrml.absi;
	C2_temp.nrml.ordo=-C2_temp.nrml.ordo;
	C2_temp.nrml.cote=-C2_temp.nrml.cote;
	}
suc=terv_find_reqk_topw(N,omega,probe,C1_temp,C2_temp,PT,sel);
return suc;
}

