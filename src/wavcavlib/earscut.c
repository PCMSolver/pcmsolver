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
#include <stdlib.h>
#include <stdio.h>
#include "cavity.h"
#include "pln_sph.h"
#include "triang.h"


int miqc_peri_quwm(int k,int N)
{int res;
res=k-1;
if(res==-1)
	res=N-1;
return res;
}


int fitw_peri_detb(int k,int N)
{int res;
res=k+1;
if(res==N)
	res=0;
return res;
}


void lumw_find_simn_hubq(parm A,parm B,parm *C,parm *D,double eps)
{C->u=(1.0-eps)*A.u+eps*B.u;
C->v=(1.0-eps)*A.v+eps*B.v;
D->u=eps*A.u+(1.0-eps)*B.u;
D->v=eps*A.v+(1.0-eps)*B.v;
}




int qolc_find_nujm_tolb(parm *P,int N,parm A,parm B,double eps)
{int i,sc,ts,res;
double marg=0.001,eps__=1.0e-11,mu=1.0e-15;
parm C,D,E,F;
lumw_find_simn_hubq(A,B,&C,&D,eps);
res=0;
for(i=0;i<N;i++)
	{sc=fitw_peri_detb(i,N);
	lumw_find_simn_hubq(P[i],P[sc],&E,&F,eps);
	
	ts=kesn_segm_lafn(C,D,E,F,marg,eps__,mu);
	if(ts==1)
		{res=1;
		break;
		}
	}
return res; 
}



int pung_test_rutp(int k,int a,int b,int c,pair_nodes PR)
{int i,res=0,f,s;
for(i=0;i<PR.nb_pairs;i++)
	{f=PR.first[i];
	s=PR.second[i];
	if((f==k)&&(s==a))
		{res=1;
		break;
		}
	if((f==k)&&(s==b))
		{res=1;
		break;
		}
	if((f==k)&&(s==c))
		{res=1;
		break;
		}
	
	if((f==a)&&(s==k))
		{res=1;
		break;
		}
	if((f==b)&&(s==k))
		{res=1;
		break;
		}
	if((f==c)&&(s==k))
		{res=1;
		break;
		}
	}
return res;
}



int cerv_segm_gusk(parm A,parm B,parm X)
{int res;
double a,b,c,eps=1.0e-6,diff;
a=pufv_dist_mekq(A,X);
b=pufv_dist_mekq(B,X);
c=pufv_dist_mekq(A,B);
diff=fabs(a+b-c);
res=0;
if(diff<eps)
	res=1;
return res;
}



void vusf_clou_kegr(parm *P,int N,double *XMIN,double *XMAX,double *YMIN,double *YMAX)
{int i;
double x,y,xmi,xma,ymi,yma;
xmi=LARGE_NUMBER;	xma=-LARGE_NUMBER;
ymi=LARGE_NUMBER;	yma=-LARGE_NUMBER;
for(i=0;i<N;i++)
	{x=P[i].u;
	y=P[i].v;
	if(x<xmi)	xmi=x;
	if(y<ymi)	ymi=y;
	if(x>xma)	xma=x;
	if(y>yma)	yma=y;
	}
*XMIN=xmi;		*XMAX=xma;
*YMIN=ymi;		*YMAX=yma;
}



int tewz_test_sowh(parm A,parm B,parm C,
double xmi,double xma,double ymi,double yma,parm P)
{int res,ts;
if(P.u<=xmi)	return 0;
if(P.u>=xma)	return 0;
if(P.v<=ymi)	return 0;
if(P.v>=yma)	return 0;

ts=cerv_segm_gusk(A,B,P);
if(ts==1)	return 1;
ts=cerv_segm_gusk(B,C,P);
if(ts==1)	return 1;
ts=cerv_segm_gusk(A,C,P);
if(ts==1)	return 1;

res=1;
ts=kujw_test_lifn(B,A,P);
if(ts==1)
	res=0;
if(res==1)
	{ts=kujw_test_lifn(C,B,P);
	if(ts==1)
		res=0;
	}
if(res==1)
	{ts=kujw_test_lifn(A,C,P);
	if(ts==1)
		res=0;
	}
return res;
}



int visq_find_jedf_fedk(int *seq,int n,parm *P,int a,int b,int c,pair_nodes PR)
{int res,i,k,ts,tr;
double xmi,xma,ymi,yma,ext=1.0e-5;
parm *temp;

temp=(parm *)malloc(3*sizeof(parm));
temp[0].u=P[a].u;		temp[0].v=P[a].v;
temp[1].u=P[b].u;		temp[1].v=P[b].v;
temp[2].u=P[c].u;		temp[2].v=P[c].v;
vusf_clou_kegr(temp,3,&xmi,&xma,&ymi,&yma);
xmi=xmi-ext;	xma=xma+ext;
ymi=ymi-ext;	yma=yma+ext;
free(temp);

res=0;
for(k=0;k<n;k++)
	{i=seq[k];
	if((i!=a)&&(i!=b)&&(i!=c))
		{tr=pung_test_rutp(i,a,b,c,PR);
		if(tr==0)
			{ts=tewz_test_sowh(P[a],P[b],P[c],xmi,xma,ymi,yma,P[i]);
			if(ts==1)
				{res=1;
				break;
				}
			}
		}
	}
return res;
}


int viqc_find_tuwm(parm Ni,parm Nj,parm Nk)
{int res;
double det;
vect R1,R2;
R1.u=Nj.u-Ni.u;   R1.v=Nj.v-Ni.v;
R2.u=Nk.u-Ni.u;   R2.v=Nk.v-Ni.v;
det=kelr_dete_lusf(R1,R2);
if(det>0.0)
	res=+1;
else
	res=-1;
return res;
}



void kuwz_circ_kilv(parm Ni,parm Nj,parm Nk,double *r_insc,double *r_circum)
{int ort;
double x,y,xij,yij,xjk,yjk,num,den,xin,yin,xi,yi,xj,yj,xk,yk;
double si, sj, sk, O,Det;
ort=viqc_find_tuwm(Ni,Nj,Nk);
if(ort==1)
	{xi=Ni.u;  yi=Ni.v;
	xj=Nj.u;  yj=Nj.v;
	xk=Nk.u;  yk=Nk.v;
	xij=0.5*(xi+xj); yij=0.5*(yi+yj);
	xjk=0.5*(xj+xk); yjk=0.5*(yj+yk);
	num = (xij-xjk)*(xj-xi) + (yij-yjk)*(yj-yi);
	den = (xj -xi) *(yk-yj) - (xk -xj) *(yj-yi);
	if(den>0)
		{x=xjk + num/den*(yk-yj);
		y=yjk - num/den*(xk-xj);
		*r_circum= sqrt( (xi-x)*(xi-x) + (yi-y)*(yi-y) );
		}
	si=pufv_dist_mekq(Ni,Nj);
	sj=pufv_dist_mekq(Nj,Nk);
	sk=pufv_dist_mekq(Nk,Ni);
	O =si+sj+sk;
	Det = xi*(yj-yk) - xj*(yi-yk) + xk*(yi-yj);
	xin = ( xi*si + xj*sj + xk*sk ) / O;
	yin = ( yi*si + yj*sj + yk*sk ) / O;
	*r_insc=Det / O;
	}
else
	kuwz_circ_kilv(Nj,Ni,Nk,r_insc,r_circum);
}



int cuql_circ_mogn(parm Ni,parm Nj,parm Nk,double *r_insc,double *r_circum)
{int suc,ort1,ort2;
ort1=viqc_find_tuwm(Ni,Nj,Nk);
ort2=viqc_find_tuwm(Nj,Ni,Nk);
if(ort1==ort2)
	suc=FAILURE;
else
	{kuwz_circ_kilv(Ni,Nj,Nk,r_insc,r_circum);
	suc=SUCCESS;
	}
return suc;
}



double hikq_aspe_qurk(parm A,parm B,parm C)
{int suc;
double r,R,res,coeff=2.0;
suc=cuql_circ_mogn(A,B,C,&r,&R);
if(suc==FAILURE)
	return LARGE_NUMBER;
res=R/(r*coeff);
return res;
}



int viqr_find_mavr_lohs(parm *P,int N,int *seq,int n,manif_ro *msh,pair_nodes PR,double eps)
{int i,k,q,pr,nx,nel,*temp,ts,suc;
double alpha,sm=0.0,rho;
suc=SUCCESS;

q=-1;
sm=10000.0;
for(i=0;i<n;i++)
	{pr=miqc_peri_quwm(i,n);
	nx=fitw_peri_detb(i,n);
	alpha=lomn_inte_cubq(P[seq[pr]],P[seq[i]],P[seq[nx]]);
	if(alpha<MY_PI)
		{rho=hikq_aspe_qurk(P[seq[pr]],P[seq[i]],P[seq[nx]]);
		fprintf(tmpout,"i=%d  rho=%f  sm=%f\n",i,rho,sm);
		if(rho<sm)
			{ts=qolc_find_nujm_tolb(P,N,P[seq[pr]],P[seq[nx]],eps);
			if(ts==0)
				{ts=visq_find_jedf_fedk(seq,n,P,seq[pr],seq[i],seq[nx],PR);
				if(ts==0)
					{sm=rho;
					q=i;
					}
				}
			}
		}
	}
if(q==-1)
	{suc=FAILURE;
	}
else
	{
	pr=miqc_peri_quwm(q,n);
	nx=fitw_peri_detb(q,n);
	nel=msh->e_grs;
	msh->entity[nel].frvrt=seq[pr];
	msh->entity[nel].scvrt=seq[q];
	msh->entity[nel].thvrt=seq[nx];
	msh->e_grs=nel+1;
	
	temp=(int *)malloc(n*sizeof(int));
	k=0;
	for(i=0;i<n;i++)
	if(i!=q)
		{temp[k]=seq[i];
		k++;
		}
	for(i=0;i<n-1;i++)
		seq[i]=temp[i];
	free(temp);
	}
return suc;
}



int rocl_find_welh_mecg(parm *P,int N,int *seq,int n,manif_ro *msh,pair_nodes PR,double eps)
{int i,k,q,pr,nx,nel,*temp,ts,suc;
double alpha;
suc=SUCCESS;

q=-1;
for(i=0;i<n;i++)
	{pr=miqc_peri_quwm(i,n);
	nx=fitw_peri_detb(i,n);
	alpha=lomn_inte_cubq(P[seq[pr]],P[seq[i]],P[seq[nx]]);
	if(alpha<MY_PI)
		{ts=qolc_find_nujm_tolb(P,N,P[seq[pr]],P[seq[nx]],eps);
		if(ts==0)
			{ts=visq_find_jedf_fedk(seq,n,P,seq[pr],seq[i],seq[nx],PR);
			if(ts==0)
				{q=i;
				break;
				}
			}
		}
	}
if(q==-1)
	suc=FAILURE;
else
	{
	pr=miqc_peri_quwm(q,n);
	nx=fitw_peri_detb(q,n);
	nel=msh->e_grs;
	msh->entity[nel].frvrt=seq[pr];
	msh->entity[nel].scvrt=seq[q];
	msh->entity[nel].thvrt=seq[nx];
	msh->e_grs=nel+1;
	
	temp=(int *)malloc(n*sizeof(int));
	k=0;
	for(i=0;i<n;i++)
	if(i!=q)
		{temp[k]=seq[i];
		k++;
		}
	for(i=0;i<n-1;i++)
		seq[i]=temp[i];
	free(temp);
	}
return suc;
}



int ear_cut_once(parm *P,int N,int *seq,
int n,manif_ro *msh,pair_nodes PR,int *forc_term)
{int suc,SUC;
double eps=1.0e-3;

*forc_term=0;
SUC=SUCCESS;
while(1)
	{if(0)
		suc=viqr_find_mavr_lohs(P,N,seq,n,msh,PR,eps);
	else
		suc=rocl_find_welh_mecg(P,N,seq,n,msh,PR,eps);
	if(suc==SUCCESS)
		break;
	if(eps<1.0e-7)
		{fprintf(tmpout,"Unable to find an ear with eps=%e\n",eps);
		
		
		*forc_term=1;
		SUC=FAILURE;
		break;
		}
	eps=0.1*eps;
	}
return SUC;
}



int fitk_find_neqs_wusr(parm *P,int N,manif_ro *msh,pair_nodes PR,
int *forc_term)
{int i,*seq,n,nel,suc,SUC,f_trm;

*forc_term=0;
for(i=0;i<N;i++)
	{msh->knot[i].u=P[i].u;
	msh->knot[i].v=P[i].v;
	}
msh->n_grs=N;
msh->e_grs=0;
n=N;
seq=(int *)malloc(N*sizeof(int));
for(i=0;i<N;i++)	
	seq[i]=i;
SUC=SUCCESS;

for(i=0;i<N;i++)
	{if(n==3)
		{nel=msh->e_grs;
		msh->entity[nel].frvrt=seq[0];
		msh->entity[nel].scvrt=seq[1];
		msh->entity[nel].thvrt=seq[2];
		msh->e_grs=nel+1;
		break;
		}
	suc=ear_cut_once(P,N,seq,n,msh,PR,&f_trm);
	if(f_trm==1)
		{fprintf(tmpout,"force term: ear_cut_once() in fitk_find_neqs_wusr()\n");
		*forc_term=1;
		free(seq);
		return FAILURE;
		}
	if(suc==FAILURE)
		{SUC=FAILURE;
		break;
		}
	n=n-1;
	}
free(seq);
return SUC;
}



int qist_list_tops(int *S,int L,int val,int *ls)
{int nb,i;
nb=0;
for(i=0;i<L;i++)
if(S[i]==val)
	{ls[nb]=i;
	nb++;
	}
return nb;
}



int zogw_pair_zalj(int *S,int L,int *val)
{int nb,i,ts,dummy;
nb=0;
for(i=0;i<L;i++)
	{ts=gonl_arra_govj(val,nb,S[i],&dummy);
	if(ts==0)
		{val[nb]=S[i];
		nb++;
		}
	}
return nb;
}



void seqd_find_wesf(int *S,int L,pair_nodes *PR)
{int N,*val,i,j,l,k,w,*ls,m;
val=(int *)malloc(L*sizeof(int));
ls=(int *)malloc(L*sizeof(int));
N=zogw_pair_zalj(S,L,val);
k=0;
for(i=0;i<N;i++)
	{w=val[i];
	m=qist_list_tops(S,L,w,ls);
	for(j=0;j<m;j++)
	for(l=0;l<j;l++)
		{PR->first[k]=ls[j];
		PR->second[k]=ls[l];
		
		k++;
		}
	}
free(val);
free(ls);
PR->nb_pairs=k;
}



void livh_find_vutj_wujn(mult_conn P,int *S,int L,manif_ro *msh)
{int n,i,n1,n2,n3;
n=P.v_grs;
for(i=0;i<n;i++)
	{msh->knot[i].u=P.vertex[i].u;
	msh->knot[i].v=P.vertex[i].v;
	}
for(i=0;i<msh->e_grs;i++)
	{n1=msh->entity[i].frvrt;	msh->entity[i].frvrt=S[n1];
	n2=msh->entity[i].scvrt;		msh->entity[i].scvrt=S[n2];
	n3=msh->entity[i].thvrt;		msh->entity[i].thvrt=S[n3];
	}
msh->n_grs=n;
}



int rofd_find_pumz_jetl(mult_conn P,manif_ro *msh,int *forc_term)
{int n,nin,*S,nnd,L,i,w,suc,f_trm;
parm *temp;
pair_nodes PR;

*forc_term=0;
n=P.v_grs;
nin=P.nb_inner_polygons;
if(nin==0)
	{temp=(parm *)malloc(n*sizeof(parm));
	PR.nb_pairs=0;
	for(i=0;i<n;i++)
		{temp[i].u=P.vertex[i].u;	
		temp[i].v=P.vertex[i].v;	
		}
	
	suc=fitk_find_neqs_wusr(temp,n,msh,PR,&f_trm);
	free(temp);
	if(f_trm==1)
		{*forc_term=1;
		fprintf(tmpout,"force term: fitk_find_neqs_wusr() in rofd_find_pumz_jetl()\n");
		return FAILURE;
		}
	}
else
	{nnd=n+2*nin;
	PR.first=(int *)malloc(n*sizeof(int));
	PR.second=(int *)malloc(n*sizeof(int));
	S=(int *)malloc(nnd*sizeof(int));
	L=qekv_gene_taqc(P,S);
	seqd_find_wesf(S,L,&PR);
	temp=(parm *)malloc(L*sizeof(parm));
	for(i=0;i<L;i++)
		{w=S[i];
		temp[i].u=P.vertex[w].u;
		temp[i].v=P.vertex[w].v;
		}
	
	suc=fitk_find_neqs_wusr(temp,nnd,msh,PR,&f_trm);
	if(f_trm==1)
		{*forc_term=1;
		fprintf(tmpout,"force term: fitk_find_neqs_wusr() in rofd_find_pumz_jetl()\n");
		free(temp);
		free(S);
		free(PR.first);
		free(PR.second);
		return FAILURE;
		}
	free(temp);
	livh_find_vutj_wujn(P,S,L,msh);
	free(S);
	free(PR.first);
	free(PR.second);
	}
return suc;
}


int vorg_find_qach_nujt(mult_conn P,manif_ro *msh,int *forc_term)
{int i,nin,suc,nnd,nel,f_trm;
double xmi,xma,ymi,yma,h_x,h_y,scl;
parm *p,G;
mult_conn temp;
manif_ro msh_temp;

*forc_term=0;
p=(parm *)malloc(P.nb_vr_outer*sizeof(parm));
for(i=0;i<P.nb_vr_outer;i++)
	cunl_find_qedf_rewn(P.vertex[i],&p[i]);
ritp_boun_niwz(p,P.nb_vr_outer,&xmi,&xma,&ymi,&yma);
G.u=0.5*(xmi+xma);	h_x=xma-xmi;
G.v=0.5*(ymi+yma);	h_y=yma-ymi;
if(h_x<h_y)
	scl=1.0/h_y;
else
	scl=1.0/h_x;
free(p);

nin=P.nb_inner_polygons;
temp.v_grs=P.v_grs;
temp.nb_vr_outer=P.nb_vr_outer;
temp.vertex=(parm *)malloc(P.v_grs*sizeof(parm));
temp.zt=(double *)malloc(P.v_grs*sizeof(double));
temp.flag=(int *)malloc(P.v_grs*sizeof(int));
temp.mapglob=(int *)malloc(P.v_grs*sizeof(int));
temp.nb_vr_inner=(int *)malloc(nin*sizeof(int));
for(i=0;i<nin;i++)
	temp.nb_vr_inner[i]=P.nb_vr_inner[i];
for(i=0;i<P.v_grs;i++)
	{temp.vertex[i].u=scl*(P.vertex[i].u-G.u);
	temp.vertex[i].v=scl*(P.vertex[i].v-G.v);
	temp.zt[i]=P.zt[i];
	temp.flag[i]=P.flag[i];
	temp.mapglob[i]=P.mapglob[i];
	}
temp.nb_inner_polygons=nin;

msh_temp.knot=(parm *)malloc(10*P.v_grs*sizeof(parm));
msh_temp.entity=(telolf *)malloc(10*P.v_grs*sizeof(telolf));
msh_temp.kt=(kt *)malloc(10*P.v_grs*sizeof(kt));
suc=rofd_find_pumz_jetl(temp,&msh_temp,&f_trm);
if(f_trm==1)
	{*forc_term=1;
	fprintf(tmpout,"force term: rofd_find_pumz_jetl() in vorg_find_qach_nujt()\n");
	free(msh_temp.knot);
	free(msh_temp.entity);
	free(msh_temp.kt);
	free(temp.vertex);
	free(temp.zt);
	free(temp.flag);
	free(temp.mapglob);
	free(temp.nb_vr_inner);
	return FAILURE;
	}
nnd=msh_temp.n_grs;
nel=msh_temp.e_grs;
for(i=0;i<nnd;i++)
	{msh->knot[i].u=(msh_temp.knot[i].u/scl)+G.u;
	msh->knot[i].v=(msh_temp.knot[i].v/scl)+G.v;
	}
for(i=0;i<nel;i++)
	{msh->entity[i].frvrt=msh_temp.entity[i].frvrt;
	msh->entity[i].scvrt=msh_temp.entity[i].scvrt;
	msh->entity[i].thvrt=msh_temp.entity[i].thvrt;
	}
msh->e_grs=nel;
msh->n_grs=nnd;

free(msh_temp.knot);
free(msh_temp.entity);
free(msh_temp.kt);
free(temp.vertex);
free(temp.zt);
free(temp.flag);
free(temp.mapglob);
free(temp.nb_vr_inner);
return suc;
}



int member_closure_triangle(parm A,parm B,parm C,parm X) 
{int ts;
double xmi,xma,ymi,yma,ext=0.1;
parm *temp;

temp=(parm *)malloc(3*sizeof(parm));
temp[0].u=A.u;		temp[0].v=A.v;
temp[1].u=B.u;		temp[1].v=B.v;
temp[2].u=C.u;		temp[2].v=C.v;
ritp_boun_niwz(temp,3,&xmi,&xma,&ymi,&yma);
xmi=xmi-ext;	xma=xma+ext;
ymi=ymi-ext;	yma=yma+ext;
free(temp);

ts=tewz_test_sowh(A,B,C,xmi,xma,ymi,yma,X);
return ts;
}
