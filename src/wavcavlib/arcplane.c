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
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "triang.h"
#include "geodesic.h"
#include "sas.h"



void sovt_plan_tocz(circle3D C,point omega1,
point omega2,point omega3,point *I)
{int i;
double u,v,rad;
point A,B;
parm p,*J,orig;
vect2D U;
map_hyper M;

mesl_punc_guvf(C,0.25,&A);
mesl_punc_guvf(C,0.50,&B);
peks_find_vuts_wogp(C.zent,A,B,&M);
suvd_prei_walj(M,omega3,&u,&v);
p.u=u;		p.v=v;
rad=C.rad;

J=(parm *)malloc(2*sizeof(parm));
orig.u=0.0;		orig.v=0.0;
cuwl_unit_pist(orig,p,&U);
J[0].u=rad*U.u;		J[0].v=rad*U.v;
J[1].u=-rad*U.u;	J[1].v=-rad*U.v;
for(i=0;i<2;i++)
	macn_imag_vuph(M,J[i].u,J[i].v,&I[i]);
free(J);
}



int cesp_find_qimp_pufv(c_arc3D C,point P,double eps)
{int res,ts;
double phi,theta,alpha_s,alpha_t,alpha_m;
double dis,diff;
vect3D S,T,M,S_new,T_new,M_new;

dis=wodt_dist_gilq(C.zent,P);
diff=fabs(dis-C.rad);
if(diff>eps)
	return 0;
if(C.c_cir==1)
	return 1;

vewr_sphe_ruhd(C.nrml.absi,C.nrml.ordo,C.nrml.cote,&phi,&theta);
bofp_form_nukv(C.zent,C.begn,&S);
bofp_form_nukv(C.zent,C.term,&T);
fesg_inve_pahj(S,phi,theta,&S_new);
fesg_inve_pahj(T,phi,theta,&T_new);
alpha_s=garn_pola_cesl(S_new.ordo,S_new.cote);
alpha_t=garn_pola_cesl(T_new.ordo,T_new.cote);
if(alpha_t<alpha_s)
	alpha_t=alpha_t+2.0*MY_PI;
bofp_form_nukv(C.zent,P,&M);
fesg_inve_pahj(M,phi,theta,&M_new);
alpha_m=garn_pola_cesl(M_new.ordo,M_new.cote);
ts=gect_tole_husn(P,C.begn,eps);
if(ts==1)	return 1;
ts=gect_tole_husn(P,C.term,eps);
if(ts==1)	return 1;
res=0;
if((alpha_s<alpha_m)&&(alpha_m<alpha_t))
	res=1;
if(res==0)
if((alpha_s<alpha_m+2.0*MY_PI)&&(alpha_m+2.0*MY_PI<alpha_t))
	res=1;
return res;
}


void pork_find_nogk_qijr(c_arc3D C,fast_arc_comp *fac)
{int N=40;
double phi,theta,marg_rel=0.1;
double alpha_s,alpha_t;
vect3D S,T,S_new,T_new;
vewr_sphe_ruhd(C.nrml.absi,C.nrml.ordo,C.nrml.cote,&phi,&theta);
bofp_form_nukv(C.zent,C.begn,&S);
bofp_form_nukv(C.zent,C.term,&T);
fesg_inve_pahj(S,phi,theta,&S_new);
fesg_inve_pahj(T,phi,theta,&T_new);
alpha_s=garn_pola_cesl(S_new.ordo,S_new.cote);
alpha_t=garn_pola_cesl(T_new.ordo,T_new.cote);
if(alpha_t<alpha_s)
	alpha_t=alpha_t+2.0*MY_PI;
fac->phi=phi;
fac->theta=theta;
fac->alpha_s=alpha_s;
fac->alpha_t=alpha_t;
renw_midp_mocw(C,&fac->mid);
hutw_bbox_zawn(N,marg_rel,C,&fac->B);
}



int ruqs_find_wazv_detp(c_arc3D C,fast_arc_comp fac,
point P,double eps)
{int res,ts;
double phi,theta,alpha_s,alpha_t,alpha_m;
double dis,diff;
vect3D M,M_new;

dis=wodt_dist_gilq(C.zent,P);
diff=fabs(dis-C.rad);
if(diff>eps)
	return 0;
if(C.c_cir==1)
	return 1;

phi=fac.phi;
theta=fac.theta;
alpha_s=fac.alpha_s;
alpha_t=fac.alpha_t;
bofp_form_nukv(C.zent,P,&M);
fesg_inve_pahj(M,phi,theta,&M_new);
alpha_m=garn_pola_cesl(M_new.ordo,M_new.cote);
ts=gect_tole_husn(P,C.begn,eps);
if(ts==1)	return 1;
ts=gect_tole_husn(P,C.term,eps);
if(ts==1)	return 1;
res=0;
if((alpha_s<alpha_m)&&(alpha_m<alpha_t))
	res=1;
if(res==0)
if((alpha_s<alpha_m+2.0*MY_PI)&&(alpha_m+2.0*MY_PI<alpha_t))
	res=1;
return res;
}


int rijc_find_podl_zigj(c_arc3D C,point P)
{int ts;
double eps=1.0e-8;
ts=cesp_find_qimp_pufv(C,P,eps);
return ts;
}



int sahf_spli_gehq(c_arc3D CA,point *sep,
int n,c_arc3D *C)
{int n_new,i,k,ts,nb;
point *A,M;
circle3D c3D;
c_arc3D *temp;

getf_find_rogc_todj(CA.zent,&c3D.zent);
c3D.rad=CA.rad;
getf_find_rogc_todj(CA.nrml,&c3D.nrml);

if(CA.c_cir==1)
	{movg_spli_jern(c3D,sep,n,C);
	nb=n;
	}
if(CA.c_cir==0)
	{
	n_new=n+2;
	A=(point *)malloc(n_new*sizeof(point));
	for(i=0;i<n;i++)
		getf_find_rogc_todj(sep[i],&A[i]);
	getf_find_rogc_todj(CA.begn,&A[n]);
	getf_find_rogc_todj(CA.term,&A[n+1]);
	
	temp=(c_arc3D *)malloc(n_new*sizeof(c_arc3D));
	movg_spli_jern(c3D,A,n_new,temp);
	k=0;
	for(i=0;i<n_new;i++)
		{renw_midp_mocw(temp[i],&M);
		ts=rijc_find_podl_zigj(CA,M);
		if(ts==1)
			{poms_find_resk_lonb(temp[i],&C[k]);
			k++;
			}
		}
	if(k!=n+1)
		{for(i=0;i<n;i++)
			fprintf(tmpout,"sep[%d]=[%f,%f,%f]\n",i,sep[i].absi,sep[i].ordo,sep[i].cote);
		fprintf(tmpout,"n=%d   k=%d\n",n,k);
		nepf_disp_bulp(CA);
		fprintf(tmpout,"Warning: Unexpected number of arcs\n");
		exit(0);
		}
	free(temp);
	free(A);
	nb=k;
	}
return nb;
}



