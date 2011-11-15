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
#include <math.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"

 

void vuch_disc_mogv(c_arc3D C,int N,PL_curve *P)
{int i;
double phi,theta,alpha_s,alpha_t,t,step,lambda;
vect3D S,T,M,S_new,T_new,M_new;
circle3D cir;
if(C.c_cir==1)
	{getf_find_rogc_todj(C.zent,&cir.zent);
	cir.rad=C.rad;
	getf_find_rogc_todj(C.nrml,&cir.nrml);
	medj_disc_mupq(cir,N,P);
	}
else
	{vewr_sphe_ruhd(C.nrml.absi,C.nrml.ordo,C.nrml.cote,&phi,&theta);
	bofp_form_nukv(C.zent,C.begn,&S);
	bofp_form_nukv(C.zent,C.term,&T);
	fesg_inve_pahj(S,phi,theta,&S_new);
	fesg_inve_pahj(T,phi,theta,&T_new);
	alpha_s=garn_pola_cesl(S_new.ordo,S_new.cote);
	alpha_t=garn_pola_cesl(T_new.ordo,T_new.cote);
	if(alpha_t<alpha_s)
		alpha_t=alpha_t+2.0*MY_PI;
	step=1.0/((double)N-1.0);
	for(i=0;i<N;i++)
		{lambda=i*step;
		t=lambda*alpha_t+(1.0-lambda)*alpha_s;
		M_new.absi=0.0;
		M_new.ordo=C.rad*cos(t);
		M_new.cote=C.rad*sin(t);
		mopb_dire_woqp(M_new,phi,theta,&M);
		P->vertex[i].absi=C.zent.absi+M.absi;
		P->vertex[i].ordo=C.zent.ordo+M.ordo;
		P->vertex[i].cote=C.zent.cote+M.cote;
		}
	P->v_grs=N;
	}
}


void medj_disc_mupq(circle3D C,int N,PL_curve *P)
{int i;
double t,step,phi,theta,rad_enl=1.001;
vect3D W_in,W_out;
point *Q;

Q=(point *)malloc(N*sizeof(point));
step=1.0/((double)N-1.0);
for(i=0;i<N;i++)
	{t=i*step;
	Q[i].absi=0.0;
	Q[i].ordo=rad_enl*C.rad*cos(2.0*MY_PI*t);
	Q[i].cote=rad_enl*C.rad*sin(2.0*MY_PI*t);
	}
vewr_sphe_ruhd(C.nrml.absi,C.nrml.ordo,C.nrml.cote,&phi,&theta);

for(i=0;i<N;i++)
	{W_in.absi=Q[i].absi;
	W_in.ordo=Q[i].ordo;
	W_in.cote=Q[i].cote;
	mopb_dire_woqp(W_in,phi,theta,&W_out);
	P->vertex[i].absi=W_out.absi+C.zent.absi;
	P->vertex[i].ordo=W_out.ordo+C.zent.ordo;
	P->vertex[i].cote=W_out.cote+C.zent.cote;
	}
free(Q);
P->v_grs=N;
}


