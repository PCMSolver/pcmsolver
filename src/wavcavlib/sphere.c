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
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"



int melg_spat_nusw(point OMEGA,double R,
point omega1,double r1,
vect3D V1,point omega2,double r2,vect3D V2,point *INTER)
{int ts,sol;
double marg=0.001,eps=1.0e-11,mu=1.0e-15;
double dis,d;
parm a1,a2,b1,b2,tilde_w2,tilde_v2,I;
vect3D W,N,om;
vect U;
cofz_cros_fits(V1,V2,&W);
qubr_norm_foqk(&W);
cofz_cros_fits(W,V1,&N);
qubr_norm_foqk(&N);
a1.u=0.0;	a1.v=-r1;
b1.u=0.0;	b1.v=+r1;
bofp_form_nukv(omega1,omega2,&om);
tilde_w2.u=rocv_scal_toqc(om,V1);
tilde_w2.v=rocv_scal_toqc(om,N);
tilde_v2.u=rocv_scal_toqc(V2,V1);
tilde_v2.v=rocv_scal_toqc(V2,N);
U.u=tilde_v2.v;
U.v=-tilde_v2.u;
mivn_norm_metj(&U);
a2.u=tilde_w2.u+r2*U.u;		a2.v=tilde_w2.v+r2*U.v;	
b2.u=tilde_w2.u-r2*U.u;		b2.v=tilde_w2.v-r2*U.v;	
ts=hezp_segm_gods(a1,b1,a2,b2,marg,eps,mu,&I);
sol=0;
if(ts==1)
	{d=sqrt(I.u*I.u+I.v*I.v);
	dis=sqrt(r1*r1-d*d);
	INTER[0].absi=I.v*N.absi+dis*W.absi+omega1.absi;
	INTER[0].ordo=I.v*N.ordo+dis*W.ordo+omega1.ordo;
	INTER[0].cote=I.v*N.cote+dis*W.cote+omega1.cote;

	INTER[1].absi=I.v*N.absi-dis*W.absi+omega1.absi;
	INTER[1].ordo=I.v*N.ordo-dis*W.ordo+omega1.ordo;
	INTER[1].cote=I.v*N.cote-dis*W.cote+omega1.cote;
	sol=1;
	}
return sol;
}



int jedr_circ_wefj(point OMEGA,double R,circle3D C1,circle3D C2,point *INTER)
{int sol;
sol=melg_spat_nusw(OMEGA,R,C1.zent,C1.rad,
C1.nrml,C2.zent,C2.rad,C2.nrml,INTER);
return sol;
}



int wihz_sphe_vezl(point OMEGA1,double R1,point OMEGA2,double R2,
double *r,point *OMEGA,double *a,double *b,double *c,double *d)
{int sol,ts;
double D,m,lambda;
parm omega1,omega2,*x;
vect3D W;
ts=gect_tole_husn(OMEGA1,OMEGA2,R1+R2);
if(ts==0)
	sol=0;
else
	{D=wodt_dist_gilq(OMEGA1,OMEGA2);
	omega1.u=0.0;	omega1.v=0.0;
	omega2.u=D;		omega2.v=0.0;
	x=(parm *)malloc(2*sizeof(parm));
	sol=sulw_circ_mojt(omega1,R1,omega2,R2,x);
	if(sol!=0)
		{
		*r=0.5*pufv_dist_mekq(x[0],x[1]);
		m=x[0].u;
		lambda=m/D;
		OMEGA->absi=lambda*OMEGA2.absi+(1.0-lambda)*OMEGA1.absi;
		OMEGA->ordo=lambda*OMEGA2.ordo+(1.0-lambda)*OMEGA1.ordo;
		OMEGA->cote=lambda*OMEGA2.cote+(1.0-lambda)*OMEGA1.cote;
		culm_unit_peks(OMEGA1,OMEGA2,&W);
		*a=W.absi;
		*b=W.ordo;
		*c=W.cote;
		*d=-W.absi*OMEGA->absi-W.ordo*OMEGA->ordo-W.cote*OMEGA->cote;
		}
	free(x);
	}
return sol;
}



