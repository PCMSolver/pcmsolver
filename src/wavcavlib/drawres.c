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
#include "splinemol.h"
#include "meshsas.h"
 

void doct_eval_setm(double r,double phi,
double theta,point omega,point *sol)
{point temp;
temp.absi=r*sin(theta)*cos(phi);
temp.ordo=r*sin(theta)*sin(phi);
temp.cote=r*cos(theta);		
sol->absi=temp.absi+omega.absi;
sol->ordo=temp.ordo+omega.ordo;
sol->cote=temp.cote+omega.cote;
}



void jirs_unit_hocd(int N,int M,int m,manif_ro *unt)
{int i,j,k_nd,k_el,n1,n2,n3;
double step_u,step_v,u,v;
double u_cur,u_nxt,v_cur,v_nxt;

step_u=1.0/((double)N*(double)m);
step_v=1.0/(double)M;
k_nd=0;		k_el=0;
for(j=0;j<=M;j++)
	{v=(double)j*step_v;
	for(i=0;i<N*m;i++)
		{u_cur=(double)i*step_u;
		u_nxt=u_cur+step_u;
		n1=k_nd;  n2=k_nd+1;  n3=k_nd+2;
		unt->knot[k_nd].u=u_cur; unt->knot[k_nd].v=v; k_nd++;
		unt->knot[k_nd].u=u_nxt; unt->knot[k_nd].v=v; k_nd++;
		unt->knot[k_nd].u=u_cur; unt->knot[k_nd].v=v; k_nd++;
		unt->entity[k_el].frvrt=n1;
		unt->entity[k_el].scvrt=n2;
		unt->entity[k_el].thvrt=n3;
		k_el++;
		}
	}

step_u=1.0/(double)N;
step_v=1.0/((double)M*(double)m);
for(j=0;j<=N;j++)
	{u=(double)j*step_u;
	for(i=0;i<M*m;i++)
		{v_cur=(double)i*step_v;
		v_nxt=v_cur+step_v;
		n1=k_nd;  n2=k_nd+1;  n3=k_nd+2;
		unt->knot[k_nd].u=u; unt->knot[k_nd].v=v_cur; k_nd++;
		unt->knot[k_nd].u=u; unt->knot[k_nd].v=v_nxt; k_nd++;
		unt->knot[k_nd].u=u; unt->knot[k_nd].v=v_cur; k_nd++;
		unt->entity[k_el].frvrt=n1;
		unt->entity[k_el].scvrt=n2;
		unt->entity[k_el].thvrt=n3;
		k_el++;
		}
	}
unt->n_grs  =k_nd;
unt->e_grs=k_el;
}



void mowb_mesh_wefs(point omega,double r,
int N,int M,int m,manif_tl *msh)
{int nnd,nel,i;
double u,v,phi,theta;
manif_ro unt;
point temp;
nel=(M+1)*N*m+(N+1)*M*m;
nnd=3*nel;
unt.knot=(parm *)malloc(nnd*sizeof(parm));
unt.entity=(telolf *)malloc(nel*sizeof(telolf));
jirs_unit_hocd(N,M,m,&unt);
nnd=unt.n_grs;
nel=unt.e_grs;
for(i=0;i<nnd;i++)
	{u=unt.knot[i].u;
	v=unt.knot[i].v;
	phi=2.0*MY_PI*u;
	theta=MY_PI*v;
	doct_eval_setm(r,phi,theta,omega,&temp);
	getf_find_rogc_todj(temp,&msh->knot[i]);
	}
for(i=0;i<nel;i++)
	{msh->entity[i].frvrt=unt.entity[i].frvrt;
	msh->entity[i].scvrt=unt.entity[i].scvrt;
	msh->entity[i].thvrt=unt.entity[i].thvrt;
	}
msh->n_grs=nnd;
msh->e_grs=nel;
free(unt.knot);
free(unt.entity);
}


void qozd_find_vijw_wont(parm A,parm B,double t,parm *sol)
{sol->u=(1.0-t)*A.u+t*B.u;
sol->v=(1.0-t)*A.v+t*B.v;
}


void korp_mapp_cafv(parm a,parm b,parm c,parm d,
double u,double v,parm *sol)
{double F0u,F0v,F1u,F1v,A,B,C;
parm alpha_u,gamma_u,delta_v,beta_v,alpha_0,alpha_1,gamma_0,gamma_1;

qozd_find_vijw_wont(a,b,u,&alpha_u);
qozd_find_vijw_wont(d,c,u,&gamma_u);
qozd_find_vijw_wont(a,d,v,&delta_v);
qozd_find_vijw_wont(b,c,v,&beta_v);

qozd_find_vijw_wont(a,b,0.0,&alpha_0);
qozd_find_vijw_wont(a,b,1.0,&alpha_1);
qozd_find_vijw_wont(d,c,0.0,&gamma_0);
qozd_find_vijw_wont(d,c,1.0,&gamma_1);
F1u=u;
F0u=1.0-F1u;
F1v=v;
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

  

void cevg_simp_cafv(int N,parm A,parm B,parm C,manif_ro *msh)
{int nnd,nel,i,NND,NEL,j,n1,n2,n3;
parm G,M1,M2,M3,a,b,c,d,temp;
manif_ro unt;
nnd=N*N;
nel=2*(N-1)*(N-1);
unt.knot=(parm *)malloc(nnd*sizeof(parm));
unt.entity=(telolf *)malloc(nel*sizeof(telolf));
lenq_simp_socg(N,N,&unt);

M1.u=0.5*(A.u+B.u);
M1.v=0.5*(A.v+B.v);

M2.u=0.5*(B.u+C.u);
M2.v=0.5*(B.v+C.v);

M3.u=0.5*(A.u+C.u);
M3.v=0.5*(A.v+C.v);

G.u=(A.u+B.u+C.u)/3.0;
G.v=(A.v+B.v+C.v)/3.0;

NND=0;	NEL=0;
for(i=0;i<3;i++)
	{if(i==0)
		{cunl_find_qedf_rewn(A,&a);
		cunl_find_qedf_rewn(M1,&b);
		cunl_find_qedf_rewn(G,&c);
		cunl_find_qedf_rewn(M3,&d);
		}
	if(i==1)
		{cunl_find_qedf_rewn(M1,&a);
		cunl_find_qedf_rewn(B,&b);
		cunl_find_qedf_rewn(M2,&c);
		cunl_find_qedf_rewn(G,&d);
		}
	if(i==2)
		{cunl_find_qedf_rewn(C,&a);
		cunl_find_qedf_rewn(M3,&b);
		cunl_find_qedf_rewn(G,&c);
		cunl_find_qedf_rewn(M2,&d);
		}
	for(j=0;j<nel;j++)
		{n1=unt.entity[j].frvrt;
		n2=unt.entity[j].scvrt;
		n3=unt.entity[j].thvrt;
		msh->entity[NEL].frvrt=NND+n1;
		msh->entity[NEL].scvrt=NND+n2;
		msh->entity[NEL].thvrt=NND+n3;
		NEL++;
		}
	for(j=0;j<nnd;j++)
		{korp_mapp_cafv(a,b,c,d,unt.knot[j].u,unt.knot[j].v,&temp);
		cunl_find_qedf_rewn(temp,&msh->knot[NND]);
		NND++;
		}
	}

msh->n_grs=NND;
msh->e_grs=NEL;
free(unt.knot);
free(unt.entity);
}



void guqw_simp_howc(int N,parm A,parm B,parm C,
manif_ro *msh,int *cr_a,int *cr_b,int *cr_c)
{int nnd,nel,i,NND,NEL,j,n1,n2,n3;
parm G,M1,M2,M3,a,b,c,d,temp;
manif_ro unt;
nnd=N*N;
nel=2*(N-1)*(N-1);
unt.knot=(parm *)malloc(nnd*sizeof(parm));
unt.entity=(telolf *)malloc(nel*sizeof(telolf));
lenq_simp_socg(N,N,&unt);

M1.u=0.5*(A.u+B.u);
M1.v=0.5*(A.v+B.v);

M2.u=0.5*(B.u+C.u);
M2.v=0.5*(B.v+C.v);

M3.u=0.5*(A.u+C.u);
M3.v=0.5*(A.v+C.v);

G.u=(A.u+B.u+C.u)/3.0;
G.v=(A.v+B.v+C.v)/3.0;

NND=0;	NEL=0;
for(i=0;i<3;i++)
	{if(i==0)
		{cunl_find_qedf_rewn(A,&a);
		cunl_find_qedf_rewn(M1,&b);
		cunl_find_qedf_rewn(G,&c);
		cunl_find_qedf_rewn(M3,&d);
		}
	if(i==1)
		{cunl_find_qedf_rewn(M1,&a);
		cunl_find_qedf_rewn(B,&b);
		cunl_find_qedf_rewn(M2,&c);
		cunl_find_qedf_rewn(G,&d);
		}
	if(i==2)
		{cunl_find_qedf_rewn(C,&a);
		cunl_find_qedf_rewn(M3,&b);
		cunl_find_qedf_rewn(G,&c);
		cunl_find_qedf_rewn(M2,&d);
		}
	for(j=0;j<nel;j++)
		{n1=unt.entity[j].frvrt;
		n2=unt.entity[j].scvrt;
		n3=unt.entity[j].thvrt;
		msh->entity[NEL].frvrt=NND+n1;
		msh->entity[NEL].scvrt=NND+n2;
		msh->entity[NEL].thvrt=NND+n3;
		NEL++;
		}
	for(j=0;j<nnd;j++)
		{korp_mapp_cafv(a,b,c,d,unt.knot[j].u,unt.knot[j].v,&temp);
		cunl_find_qedf_rewn(temp,&msh->knot[NND]);
		if((i==0)&&(j==0))
			*cr_a=NND;
		if((i==1)&&(j==N-1))
			*cr_b=NND;
		if((i==2)&&(j==0))
			*cr_c=NND;
		NND++;
		}
	}

msh->n_grs=NND;
msh->e_grs=NEL;
free(unt.knot);
free(unt.entity);
}

  

void jatw_unit_hukl(int N,manif_ro *msh)
{parm A,B,C;
A.u=0.0;	A.v=0.0;
B.u=1.0;	B.v=0.0;
C.u=0.0;	C.v=1.0;
cevg_simp_cafv(N,A,B,C,msh);
}


int wucl_find_henl_lorh(trmsrf surf)
{
return 0;
}
 

void jegn_find_wikc_poqc(manif_tl msh,int proc)
{int nnd,i;
FILE *fp;
if(proc==0)	fp=fopen("MIOTY/grd_nodes.dat","w");
if(proc==1)	fp=fopen("MIOTY/trm_nodes.dat","w");
nnd=msh.n_grs;
for(i=0;i<nnd;i++)
	fprintf(fp,"%f  %f   %f\n",msh.knot[i].absi,msh.knot[i].ordo,msh.knot[i].cote);
fclose(fp);
}


void rinf_find_qijc_wuvr(manif_tl msh,int proc)
{int nel,i,n1,n2,n3;
FILE *fp;
if(proc==0)	fp=fopen("MIOTY/grd_elements.dat","w");
if(proc==1)	fp=fopen("MIOTY/trm_elements.dat","w");
nel=msh.e_grs;
for(i=0;i<nel;i++)
	{n1=msh.entity[i].frvrt;
	n2=msh.entity[i].scvrt;
	n3=msh.entity[i].thvrt;
	fprintf(fp,"%d  %d  %d\n",n1+1,n2+1,n3+1);
	}
fclose(fp);
}


void vumf_find_pobs_gafp(manif_tl msh,int proc)
{int nel,i;
FILE *fp;
if(proc==0)	fp=fopen("MIOTY/grd_tcolor.dat","w");
if(proc==1)	fp=fopen("MIOTY/trm_tcolor.dat","w");
nel=msh.e_grs;
for(i=0;i<nel;i++)
	{if(proc==0)	
		fprintf(fp,"0.0  0.0  0.0\n");
	if(proc==1)	
		fprintf(fp,"0.0  1.0  0.0\n");
	}
fclose(fp);
}

 
void josh_expo_qazg(point omega,double rad)
{int N=15,M=15,m=5,nnd,nel;
manif_tl msh;
nel=(M+1)*N*m+(N+1)*M*m;
nnd=3*nel;
msh.knot=(point *)malloc(nnd*sizeof(point));
msh.entity=(telolf *)malloc(nel*sizeof(telolf));
mowb_mesh_wefs(omega,rad,N,M,m,&msh);
jegn_find_wikc_poqc(msh,0);
rinf_find_qijc_wuvr(msh,0);
vumf_find_pobs_gafp(msh,0);

free(msh.entity);
free(msh.knot);
}







