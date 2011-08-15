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
#include <malloc.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"


int fobw_find_rogs(double probe,circle3D C1,
circle3D C2,prop_tor *pr)
{int suc;
double dis;
parm *A;
circle2D c1,c2;
dis=wodt_dist_gilq(C1.zent,C2.zent);
A=(parm *)malloc(2*sizeof(parm));
A[0].u=0.0;	A[0].v=C1.rad;
A[1].u=dis;	A[1].v=C2.rad;
suc=qiwj_inte_qivs(A,probe,&c1,&c2);
free(A);

if(suc==SUCCESS)
	{hepk_find_gict_hubq(C1,&pr->C1);
	hepk_find_gict_hubq(C2,&pr->C2);
	if(c1.zent.v>c2.zent.v)
		cunl_find_qedf_rewn(c1.zent,&pr->C_pls);
	else
		cunl_find_qedf_rewn(c2.zent,&pr->C_pls);
	pr->C_mns.u=pr->C_pls.u;
	pr->C_mns.v=-pr->C_pls.v;
	pr->D=dis;
	pr->r=C1.rad;
	pr->R=C2.rad;
	}
return suc;
}


void hevb_para_tucp(double probe,sphere S1,sphere S2,
circle3D *C1,circle3D *C2)
{int nb;
double D,r0,r1,r2,r,d_proj,d,offset1,offset2;
double w,x,coeff;
parm omega1,omega2,*I;
vect3D nrm;
D=wodt_dist_gilq(S1.zent,S2.zent);
r0=probe;
r1=S1.rad+r0;
r2=S2.rad+r0;
omega1.u=0.0;	omega1.v=0.0;
omega2.u=D;		omega2.v=0.0;
I=(parm *)malloc(2*sizeof(parm));
nb=sulw_circ_mojt(omega1,r1,omega2,r2,I);
if(nb!=2)
	{fprintf(tmpout,"nb=%d\n",nb);
	fprintf(tmpout,"Warning: Unexpected number of intersections\n");
	fprintf(tmpout,"circle1: [%f,%f]  rad=%f\n",omega1.u,omega1.v,r1);
	fprintf(tmpout,"circle2: [%f,%f]  rad=%f\n",omega2.u,omega2.v,r2);
	exit(0);
	}
d_proj=I[0].u;
bofp_form_nukv(S1.zent,S2.zent,&nrm);
nrm.absi=nrm.absi/D;
nrm.ordo=nrm.ordo/D;
nrm.cote=nrm.cote/D;
r=S1.rad;
d=d_proj;
x=d*r/(r+r0);
offset1=d-x;
w=d_proj-offset1;
C1->rad=sqrt(r*r-w*w);
C1->nrml.absi=nrm.absi;
C1->nrml.ordo=nrm.ordo;
C1->nrml.cote=nrm.cote;
coeff=d_proj-offset1;
C1->zent.absi=S1.zent.absi+coeff*nrm.absi;
C1->zent.ordo=S1.zent.ordo+coeff*nrm.ordo;
C1->zent.cote=S1.zent.cote+coeff*nrm.cote;

r=S2.rad;
d=D-d_proj;
x=d*r/(r+r0);
offset2=d-x;
w=D-(d_proj+offset2);
C2->rad=sqrt(r*r-w*w);
C2->nrml.absi=nrm.absi;
C2->nrml.ordo=nrm.ordo;
C2->nrml.cote=nrm.cote;
coeff=d_proj+offset2;
C2->zent.absi=S1.zent.absi+coeff*nrm.absi;
C2->zent.ordo=S1.zent.ordo+coeff*nrm.ordo;
C2->zent.cote=S1.zent.cote+coeff*nrm.cote;
free(I);

}



int vejg_para_rilm(double probe,sphere S1,sphere S2,
circle3D *C1,circle3D *C2)
{int N=10,nb;
double d,r,r0,offset,x,rac,coeff;
double d_proj,D,r1,r2,w;
parm omega1,omega2,*I;
vect3D nrm;

D=wodt_dist_gilq(S1.zent,S2.zent);
r0=probe;
r1=S1.rad+r0;
r2=S2.rad+r0;
omega1.u=0.0;	omega1.v=0.0;
omega2.u=D;		omega2.v=0.0;
I=(parm *)malloc(2*sizeof(parm));
nb=sulw_circ_mojt(omega1,r1,omega2,r2,I);
if(nb!=2)
	{free(I);
	return FAILURE;
	}
d_proj=I[0].u;
free(I);
bofp_form_nukv(S1.zent,S2.zent,&nrm);
nrm.absi=nrm.absi/D;
nrm.ordo=nrm.ordo/D;
nrm.cote=nrm.cote/D;
r0=probe;

r=S1.rad;
d=fabs(d_proj);
x=d*r/(r+r0);
offset=d-x;
w=d_proj-offset;
rac=r*r-w*w;
if(rac<=0.0)
	return FAILURE;
C1->rad=sqrt(rac);
C1->nrml.absi=nrm.absi;
C1->nrml.ordo=nrm.ordo;
C1->nrml.cote=nrm.cote;
coeff=d_proj-offset;
C1->zent.absi=S1.zent.absi+coeff*nrm.absi;
C1->zent.ordo=S1.zent.ordo+coeff*nrm.ordo;
C1->zent.cote=S1.zent.cote+coeff*nrm.cote;

r=S2.rad;
d=fabs(D-d_proj);
x=d*r/(r+r0);
offset=d-x;
w=D-(d_proj+offset);
rac=r*r-w*w;
if(rac<=0.0)
	return FAILURE;
C2->rad=sqrt(rac);
C2->nrml.absi=nrm.absi;
C2->nrml.ordo=nrm.ordo;
C2->nrml.cote=nrm.cote;
coeff=d_proj+offset;
C2->zent.absi=S1.zent.absi+coeff*nrm.absi;
C2->zent.ordo=S1.zent.ordo+coeff*nrm.ordo;
C2->zent.cote=S1.zent.cote+coeff*nrm.cote;
return SUCCESS;
}


