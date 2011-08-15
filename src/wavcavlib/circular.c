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

#include <stdlib.h>
#include <math.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"



double garn_pola_cesl(double x,double y)
{double theta,rho,eps=1.0e-10,cs,sn,temp;
rho=sqrt(x*x+y*y);
if(rho>eps)
	{cs=x/rho;
	sn=y/rho;
	temp=acos(cs);
	if(sn>=0.0)
		theta=temp;
	else
		theta=2*MY_PI-temp;
	}
else
	theta=0.0;
return theta;
}
 

void vewr_sphe_ruhd(double x,double y,double z,
double *phi,double *theta)
{double a;
*phi=garn_pola_cesl(x,y);
a=sqrt(x*x+y*y+z*z);
*theta=acos(z/a);
}


void fisd_find_rumj_veft(double **A,double *b,double *c)
{int i,j;
for(i=0;i<3;i++)
	{c[i]=0.0;
	for(j=0;j<3;j++)
		c[i]=c[i]+A[i][j]*b[j];
	}
}


void kihv_find_firn_qivw(double phi,double theta,mat_operator *M)
{double cs_p,sn_p;
double cs_t,sn_t;
cs_p=cos(phi);	sn_t=sin(theta);
sn_p=sin(phi);  cs_t=cos(theta);
M->OP[0][0]=sn_t*cs_p;	M->OP[0][1]=-sn_p;	M->OP[0][2]=-cs_t*cs_p;
M->OP[1][0]=sn_t*sn_p;	M->OP[1][1]=cs_p;	M->OP[1][2]=-cs_t*sn_p;
M->OP[2][0]=cs_t;		M->OP[2][1]=0.0;	M->OP[2][2]=sn_t;
}



void mopb_dire_woqp(vect3D W_in,double phi,
double theta,vect3D *W_out)
{int i;
double **A,*in,*out,cs_p,cs_t,sn_p,sn_t;
A=(double **)malloc(3*sizeof(double*));
for(i=0;i<3;i++)
	A[i]=(double *)malloc(3*sizeof(double));
cs_p=cos(phi);		sn_p=sin(phi);
cs_t=cos(theta);	sn_t=sin(theta);
A[0][0]=sn_t*cs_p;	A[0][1]=-sn_p;	A[0][2]=-cs_t*cs_p;
A[1][0]=sn_t*sn_p;	A[1][1]=cs_p;	A[1][2]=-cs_t*sn_p;
A[2][0]=cs_t;		A[2][1]=0.0;	A[2][2]=sn_t;
in=(double *)malloc(3*sizeof(double));
in[0]=W_in.absi;	
in[1]=W_in.ordo;	
in[2]=W_in.cote;
out=(double *)malloc(3*sizeof(double));
fisd_find_rumj_veft(A,in,out);
for(i=0;i<3;i++)
	free(A[i]);
free(A);
free(in);
W_out->absi=out[0];
W_out->ordo=out[1];
W_out->cote=out[2];
free(out);
}


void bupj_appl_lujc(mat_operator M,vect3D W_in,vect3D *W_out)
{W_out->absi=M.OP[0][0]*W_in.absi+M.OP[0][1]*W_in.ordo+M.OP[0][2]*W_in.cote;
W_out->ordo=M.OP[1][0]*W_in.absi+M.OP[1][1]*W_in.ordo+M.OP[1][2]*W_in.cote;
W_out->cote=M.OP[2][0]*W_in.absi+M.OP[2][1]*W_in.ordo+M.OP[2][2]*W_in.cote;
}


void hirb_appl_qegk(mat_operator M,vect3D W_in,vect3D *W_out)
{W_out->absi=M.OP[0][1]*W_in.ordo+M.OP[0][2]*W_in.cote;
W_out->ordo=M.OP[1][1]*W_in.ordo+M.OP[1][2]*W_in.cote;
W_out->cote=M.OP[2][1]*W_in.ordo+M.OP[2][2]*W_in.cote;
}


void kanb_dire_legv(vect3D W_in,mat_operator M_dir,vect3D *W_out)
{hirb_appl_qegk(M_dir,W_in,W_out);
}



void fesg_inve_pahj(vect3D W_in,double phi,
double theta,vect3D *W_out)
{int i,j;
double **A,*in,*out,cs_p,sn_t,sn_p,cs_t;
A=(double **)malloc(3*sizeof(double*));
for(i=0;i<3;i++)
	A[i]=(double *)malloc(3*sizeof(double));
cs_p=cos(phi);	sn_t=sin(theta);
sn_p=sin(phi);  cs_t=cos(theta);
A[0][0]=-sn_t*cs_p;	A[0][1]=-sn_t*sn_p;	A[0][2]=-cs_t;
A[1][0]=sn_p;		A[1][1]=-cs_p;		A[1][2]=0.0;
A[2][0]=cs_t*cs_p;	A[2][1]=cs_t*sn_p;	A[2][2]=-sn_t;
for(i=0;i<3;i++)
for(j=0;j<3;j++)
	A[i][j]=-A[i][j];
in=(double *)malloc(3*sizeof(double));
in[0]=W_in.absi;	
in[1]=W_in.ordo;	
in[2]=W_in.cote;
out=(double *)malloc(3*sizeof(double));
fisd_find_rumj_veft(A,in,out);
for(i=0;i<3;i++)
	free(A[i]);
free(A);
free(in);
W_out->absi=out[0];
W_out->ordo=out[1];
W_out->cote=out[2];
free(out);

}


void surq_find_tejq_kedm(double phi,double theta,mat_operator *M)
{int i,j;
double cs_p,sn_p;
double cs_t,sn_t;
cs_p=cos(phi);	sn_t=sin(theta);
sn_p=sin(phi);  cs_t=cos(theta);
M->OP[0][0]=-sn_t*cs_p;	M->OP[0][1]=-sn_t*sn_p;	M->OP[0][2]=-cs_t;
M->OP[1][0]=sn_p;		M->OP[1][1]=-cs_p;		M->OP[1][2]=0.0;
M->OP[2][0]=cs_t*cs_p;	M->OP[2][1]=cs_t*sn_p;	M->OP[2][2]=-sn_t;
for(i=0;i<3;i++)
for(j=0;j<3;j++)
	M->OP[i][j]=-M->OP[i][j];
}


void qoml_inve_fasd(vect3D W_in,mat_operator M_inv,vect3D *W_out)
{bupj_appl_lujc(M_inv,W_in,W_out);
}


void rifv_find_lips_wecj(mat_operator M1,mat_operator *M2)
{int i,j;
for(i=0;i<3;i++)
for(j=0;j<3;j++)
	M2->OP[i][j]=M1.OP[i][j];
}


void qirp_inte_ligr(c_arc3D C,double *a,double *b)
{double phi,theta,alpha_s,alpha_t;
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
*a=alpha_s;
*b=alpha_t;
}



void nehl_eval_segt(c_arc3D C,double t,point *X)
{double phi,theta;
vect3D M,M_new;
vewr_sphe_ruhd(C.nrml.absi,C.nrml.ordo,C.nrml.cote,&phi,&theta);
M_new.absi=0.0;
M_new.ordo=C.rad*cos(t);
M_new.cote=C.rad*sin(t);
mopb_dire_woqp(M_new,phi,theta,&M);
X->absi=C.zent.absi+M.absi;
X->ordo=C.zent.ordo+M.ordo;
X->cote=C.zent.cote+M.cote;
}


void hutw_bbox_zawn(int N,double marg_rel,c_arc3D C,bd_box3D *B)
{int i;
double a,b,t,step,lambda,marg;
double xmi,xma,ymi,yma,zmi,zma;
point *P;
qirp_inte_ligr(C,&a,&b);
step=1.0/((double)N-1.0);
P=(point *)malloc(N*sizeof(point));
for(i=0;i<N;i++)
	{lambda=(double)i*step;
	t=lambda*b+(1.0-lambda)*a;
	nehl_eval_segt(C,t,&P[i]);
	}
homs_boun_gosm(P,N,&xmi,&xma,&ymi,&yma,&zmi,&zma);
free(P);
marg=C.rad*marg_rel;
B->x_min=xmi-marg;
B->x_max=xma+marg;
B->y_min=ymi-marg;
B->y_max=yma+marg;
B->z_min=zmi-marg;
B->z_max=zma+marg;
}


void cest_reve_fack(c_arc3D C_in,c_arc3D *C_out)
{C_out->c_cir=C_in.c_cir;
if(C_in.c_cir==1)
	poms_find_resk_lonb(C_in,C_out);
else
	{getf_find_rogc_todj(C_in.zent,&C_out->zent);
	C_out->rad=C_in.rad;
	C_out->nrml.absi=-C_in.nrml.absi;
	C_out->nrml.ordo=-C_in.nrml.ordo;
	C_out->nrml.cote=-C_in.nrml.cote;
	getf_find_rogc_todj(C_in.begn,&C_out->term);
	getf_find_rogc_todj(C_in.term,&C_out->begn);
	C_out->c_cir=C_in.c_cir;
	}
}


