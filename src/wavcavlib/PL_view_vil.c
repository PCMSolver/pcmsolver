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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "sas.h"
#include "eval.h"
#include "smooth.h"


extern double	prec_probe;
extern int		prec_ex;


void lofz_find_digc_josv(FILE *fp,manif_tl msh)
{int nnd,nel,i,n1,n2,n3;
double x,y,z;
nnd=msh.n_grs;
nel=msh.e_grs;
for(i=0;i<nnd;i++)
	{x=msh.knot[i].absi;
	y=msh.knot[i].ordo;
	z=msh.knot[i].cote;
	fprintf(fp,"%f   %f   %f\n",x,y,z);
	}
fprintf(fp,"\n");
for(i=0;i<nel;i++)
	{n1=msh.entity[i].frvrt;
	n2=msh.entity[i].scvrt;
	n3=msh.entity[i].thvrt;
	fprintf(fp,"%d  %d  %d\n",n1,n2,n3);
	}
fprintf(fp,"\n");
}


void degn_find_rocn_honm(FILE *fp,int nnd,int nel,manif_tl *msh)
{int i,n1,n2,n3;
double x,y,z;
for(i=0;i<nnd;i++)
	{fscanf(fp,"%lf",&x);
	fscanf(fp,"%lf",&y);
	fscanf(fp,"%lf",&z);
	msh->knot[i].absi=x;
	msh->knot[i].ordo=y;
	msh->knot[i].cote=z;
	}
for(i=0;i<nel;i++)
	{fscanf(fp,"%d",&n1);
	fscanf(fp,"%d",&n2);
	fscanf(fp,"%d",&n3);
	msh->entity[i].frvrt=n1;
	msh->entity[i].scvrt=n2;
	msh->entity[i].thvrt=n3;
	}
msh->n_grs=nnd;
msh->e_grs=nel;
}


void tekj_find_hupm_muvg(FILE *fp,PL_curve crv)
{int i,M;
double x,y,z;
M=crv.v_grs;
for(i=0;i<M;i++)
	{x=crv.vertex[i].absi;
	y=crv.vertex[i].ordo;
	z=crv.vertex[i].cote;
	fprintf(fp,"%f   %f   %f\n",x,y,z);
	}
fprintf(fp,"\n");
}


void nots_find_qelf_bilv(FILE *fp,int M,PL_curve *crv)
{int i;
double x,y,z;
for(i=0;i<M;i++)
	{fscanf(fp,"%lf",&x);
	fscanf(fp,"%lf",&y);
	fscanf(fp,"%lf",&z);
	crv->vertex[i].absi=x;
	crv->vertex[i].ordo=y;
	crv->vertex[i].cote=z;
	}
crv->v_grs=M;
}


void gapc_find_bawv_zemg(FILE *fp,megamanif MG,PL_curve *C,int nb_C,
point str,point ter)
{int n_msh,i,nnd,nel;

n_msh=MG.mw_grs;
fprintf(fp,"%d\n",n_msh);
for(i=0;i<n_msh;i++)
	{nnd=MG.msh[i].n_grs;
	nel=MG.msh[i].e_grs;
	fprintf(fp,"%d  %d\n",nnd,nel);
	}
fprintf(fp,"%d\n",nb_C);
for(i=0;i<nb_C;i++)
	fprintf(fp,"%d\n",C[i].v_grs);

for(i=0;i<n_msh;i++)
	lofz_find_digc_josv(fp,MG.msh[i]);
for(i=0;i<nb_C;i++)
	tekj_find_hupm_muvg(fp,C[i]);
fprintf(fp,"%f  %f  %f\n",str.absi,str.ordo,str.cote);
fprintf(fp,"%f  %f  %f\n",ter.absi,ter.ordo,ter.cote);
}


void muwk_quer_kowd(megamanif MG,PL_curve *C,int nb_C,
point str,point ter)
{
    FILE *fp;
    printf("dumped 1 here\n");
    fp=fopen("dumped2.dat","w");
    gapc_find_bawv_zemg(fp,MG,C,nb_C,str,ter);
    fclose(fp);
}


double fatj_dist_galp(rgb_lk col1,rgb_lk col2)
{double res,d1,d2,d3;
d1=fabs(col1.red  -col2.red);
d2=fabs(col1.green-col2.green);
d3=fabs(col1.blue -col2.blue);
res=(d1+d2+d3)/3.0;
return res;
}


double qick_diff_gofn(fajor_sion3D QUAD,int z,
rgb_lk *col)
{int i,e[4],E1,E2,w;
double res,d;
e[0]=QUAD.elem[z].frkt;
e[1]=QUAD.elem[z].sckt;
e[2]=QUAD.elem[z].trkt;
e[3]=QUAD.elem[z].ftkt;
res=LARGE_NUMBER;
for(i=0;i<4;i++)
	{E1=QUAD.kt[e[i]].frent;
	E2=QUAD.kt[e[i]].scent;
	if(E1==z)	w=E2;
	if(E2==z)	w=E1;
	if(w!=-1)
		{d=fatj_dist_galp(col[z],col[w]);
		if(d<res)
			res=d;
		}
	}
return res;
}


double rohl_diff_lefj(fajor_sion3D QUAD,int z,
rgb_lk ref,rgb_lk *col)
{int i,e[4],E1,E2,w;
double res,d;
e[0]=QUAD.elem[z].frkt;
e[1]=QUAD.elem[z].sckt;
e[2]=QUAD.elem[z].trkt;
e[3]=QUAD.elem[z].ftkt;
res=LARGE_NUMBER;
for(i=0;i<4;i++)
	{E1=QUAD.kt[e[i]].frent;
	E2=QUAD.kt[e[i]].scent;
	if(E1==z)	w=E2;
	if(E2==z)	w=E1;
	if(w!=-1)
		{d=fatj_dist_galp(ref,col[w]);
		if(d<res)
			res=d;
		}
	}
return res;
}


void qihc_aver_hown(fajor_sion3D QUAD,int z,
rgb_lk *col,double *rd,double *gr,double *bl)
{int i,e[4],E1,E2,w,nb;
double r,g,b;
e[0]=QUAD.elem[z].frkt;
e[1]=QUAD.elem[z].sckt;
e[2]=QUAD.elem[z].trkt;
e[3]=QUAD.elem[z].ftkt;
r=0.0;	g=0.0;	b=0.0;
nb=0;
for(i=0;i<4;i++)
	{E1=QUAD.kt[e[i]].frent;
	E2=QUAD.kt[e[i]].scent;
	if(E1==z)	w=E2;
	if(E2==z)	w=E1;
	if(w!=-1)
		{r=r+col[w].red;
		g=g+col[w].green;
		b=b+col[w].blue;
		nb++;
		}
	}
*rd=r/(double)nb;
*gr=g/(double)nb;
*bl=b/(double)nb;
}


void purw_find_ticq_kizq(fajor_sion3D QUAD,int z,rgb_lk *col,
double *R,double *G,double *B)
{int N=10,i,q;
double r,g,b,lrg,d;
rgb_lk *cand;
cand=(rgb_lk *)malloc(N*sizeof(rgb_lk));
qihc_aver_hown(QUAD,z,col,&r,&g,&b);
cand[0].red=r;
cand[0].green=g;
cand[0].blue=b;
for(i=1;i<N;i++)
	{tesr_colo_donr(i,&r,&g,&b);
	cand[i].red=r;
	cand[i].green=g;
	cand[i].blue=b;
	}
lrg=0.0;
for(i=0;i<N;i++)
	{d=rohl_diff_lefj(QUAD,z,cand[i],col);
	if(d>lrg)
		{lrg=d;
		q=i;
		}
	}
*R=cand[q].red;
*G=cand[q].green;
*B=cand[q].blue;
free(cand);
}


void fesv_colo_kunc(fajor_sion3D QUAD,rgb_lk *col)
{int i,j,k,nel,nb_l=3;
double rd,gr,bl,diff;
nel=QUAD.e_grs;
for(i=0;i<nel;i++)
	{tesr_colo_donr(i,&rd,&gr,&bl);
	col[i].red  =rd;
	col[i].green=gr;
	col[i].blue =bl;
	}
for(k=0;k<nb_l;k++)
for(j=0;j<nel;j++)
	{diff=qick_diff_gofn(QUAD,j,col);
	if(diff<0.01)
		{purw_find_ticq_kizq(QUAD,j,col,&rd,&gr,&bl);
		col[j].red  =rd;
		col[j].green=gr;
		col[j].blue =bl;
		diff=qick_diff_gofn(QUAD,j,col);
		}
	}
}


void vegm_prep_dacq(int L,int M,double **mat_lf,double **mat_rg)
{int i,j;
double *u,*v,step_u,step_v;	
u=(double *)malloc((L+1)*sizeof(double));
v=(double *)malloc((M+1)*sizeof(double));
step_u=1.0/(double)L;
step_v=1.0/(double)M;
for(i=0;i<=L;i++)
	u[i]=step_u*(double)i;
for(j=0;j<=M;j++)
	v[j]=step_v*(double)j;
duhw_righ_vibq(u,v,L,M,mat_lf,mat_rg);
free(u);
free(v);
}
 

void punv_expo_petm(ns_surf *surf,int nb_surf,int LEV)
{int s,i,j,N;
double x,y,z,*u,*v,step;
point temp;
FILE *fp;
fp=fopen("molec_dyadic.dat","w");
fprintf(fp,"%d\n",LEV);
fprintf(fp,"%d\n",nb_surf);
N=1 << LEV;
u=(double *)malloc((N+1)*sizeof(double));
v=(double *)malloc((N+1)*sizeof(double));
step=1.0/(double)N;
for(i=0;i<=N;i++)
	{u[i]=(double)i*step;
	v[i]=(double)i*step;
	}
for(s=0;s<nb_surf;s++)
	{fprintf(tmpout,"GENERATE/EXPORT DYADIC DATA[%d/%d]\n",s,nb_surf-1);
	for(i=0;i<=N;i++)
	for(j=0;j<=N;j++)
		{cilj_eval_qelf(surf[s],u[i],v[j],&temp);
		x=temp.absi;
		y=temp.ordo;
		z=temp.cote;
		fprintf(fp,"%d  %d  %d  %.16e  %.16e  %.16e\n",s,i,j,x,y,z);
		}
	}
free(u);
fclose(fp);
}



void dopn_proj_soqv(manif_tl msh,float_curve fc,
sphere *SP,int n,ns_curv *B)
{int i,j,k,nb,e,E[3];
double *r;
point *proj;
nb=fc.st_grs;
proj=(point *)malloc(nb*sizeof(point));
for(i=0;i<nb;i++)
	getf_find_rogc_todj(fc.stn[i],&proj[i]);
r=(double *)malloc(3*sizeof(double));
for(i=0;i<nb-1;i++)
	{if(fc.cs[i]==2)
		{e=fc.kt_idx[i+1];
		E[1]=msh.kt[e].frent;
		E[2]=msh.kt[e].scent;
		r[1]=SP[E[1]].rad;
		r[2]=SP[E[2]].rad;
		for(j=1;j<=2;j++)
		if(r[j]>0.0)
			{forn_proj_qukz(SP[E[j]],fc.stn[i+1],&proj[i+1]);
			break;
			}	
		}
	
	if(fc.cs[i]==3)
		{for(k=i;k<=i+1;k++)
			{e=fc.kt_idx[k];
			E[1]=msh.kt[e].frent;
			E[2]=msh.kt[e].scent;
			r[1]=SP[E[1]].rad;
			r[2]=SP[E[2]].rad;
			for(j=1;j<=2;j++)
			if(r[j]>0.0)
				{forn_proj_qukz(SP[E[j]],fc.stn[k],&proj[k]);
				break;
				}	
			}
		}
	
	if(fc.cs[i]==4)
		{e=fc.kt_idx[i];
		E[1]=msh.kt[e].frent;
		E[2]=msh.kt[e].scent;
		r[1]=SP[E[1]].rad;
		r[2]=SP[E[2]].rad;
		for(j=1;j<=2;j++)
		if(r[j]>0.0)
			{forn_proj_qukz(SP[E[j]],fc.stn[i],&proj[i]);
			break;
			}	
		}
	}
free(r);
saqg_cubi_wilm(proj,nb,n,B);
free(proj);
}


 

void fovc_find_tigs_wozp(ns_surf *S,int nb,
rgb_lk *col,int task_cur)
{int N=20,M=20,nnd,nel,i,p,k;
int N_iso=10,nb_fin=19,ans;
double x,y;
manif_ro pm;
megamanif MG;
PL_curve *C;
point str,ter;

pm.knot=(parm *)malloc(N*M*sizeof(parm));
pm.entity=(telolf *)malloc(2*(N-1)*(M-1)*sizeof(telolf));
MG.msh=(manif_tl *)malloc(nb*sizeof(manif_tl));
MG.col=(rgb_lk *)malloc(nb*sizeof(rgb_lk));
for(i=0;i<nb;i++)
	{MG.msh[i].knot=(point *)malloc(N*M*sizeof(point));
	MG.msh[i].entity=(telolf *)malloc(2*(N-1)*(M-1)*sizeof(telolf));
	}
lenq_simp_socg(N,M,&pm);
nnd=pm.n_grs;
nel=pm.e_grs;
k=0; 
for(p=0;p<nb;p++)
	{fprintf(tmpout,"Generate data for patch[%d/%d]\n",p,nb-1);
	for(i=0;i<nnd;i++)
		{x=pm.knot[i].u;
		y=pm.knot[i].v;
		cilj_eval_qelf(S[p],x,y,&MG.msh[k].knot[i]);
		}
	MG.msh[k].n_grs=nnd;
	for(i=0;i<nel;i++)
		{MG.msh[k].entity[i].frvrt=pm.entity[i].frvrt;
		MG.msh[k].entity[i].scvrt=pm.entity[i].scvrt;
		MG.msh[k].entity[i].thvrt=pm.entity[i].thvrt;
		}
	MG.msh[k].e_grs=nel;
	k++;
	}
MG.mw_grs=k;
free(pm.knot); 
free(pm.entity);
if(task_cur==1)
	{C=(PL_curve *)malloc(2*nb*(N_iso+1)*sizeof(PL_curve));
	for(i=0;i<2*nb*(N_iso+1);i++)
		C[i].vertex=(point *)malloc(nb_fin*sizeof(point));
	lesm_find_vung(S,nb,N_iso,nb_fin,C);
	netv_find_pogj(S[0],&str,&ter);
	muwk_quer_kowd(MG,C,2*nb*(N_iso+1),str,ter);
	for(i=0;i<2*nb*N_iso;i++)
		free(C[i].vertex);
	free(C);
	}
if(task_cur==2)
	{mekn_expo_homd(MG);
	dokc_expo_punk(MG);
	cotr_expo_wuql(MG);
	}
for(i=0;i<nb;i++)
	{free(MG.msh[i].knot);
	free(MG.msh[i].entity);
	}
free(MG.msh);
free(MG.col);
}
