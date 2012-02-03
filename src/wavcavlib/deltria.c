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
#include "pln_sph.h"
#include "triang.h"



double lomf_angl_serb(parm A,parm B,parm C)
{double res,nm,su,sv,tu,tv,sp;
nm=sqrt((B.u-A.u)*(B.u-A.u)+(B.v-A.v)*(B.v-A.v));
su=(B.u-A.u)/nm;   sv=(B.v-A.v)/nm;
nm=sqrt((C.u-A.u)*(C.u-A.u)+(C.v-A.v)*(C.v-A.v));
tu=(C.u-A.u)/nm;   tv=(C.v-A.v)/nm;
sp=su*tu+sv*tv;
res=acos(sp);
return res;
}



int sujm_test_fujk(kt *ed,int sn,int tn,int N,int *flag)
{int i;
*flag=-1;
for(i=0;i<N;i++)
	{if((ed[i].frvrt==sn)&&(ed[i].scvrt==tn))
		{*flag=1;
		break;
		}
	if((ed[i].frvrt==tn)&&(ed[i].scvrt==sn))
		{*flag=1;
		break;
		}
	}
return i;
}



int vawc_dete_piln(manif_ro msh,kt *ed,int **T)
{int i,nnd,nel,n1,n2,n3,k,N,s,t,j;
nnd=msh.n_grs;
nel=msh.e_grs;
k=0;N=0;
for(i=0;i<nel;i++)
	{j=1;
	n1=msh.entity[i].frvrt;
	n2=msh.entity[i].scvrt;
	n3=msh.entity[i].thvrt;
	t=sujm_test_fujk(ed,n1,n2,N,&s);
	if(s==-1)
		{N++;
		ed[k].frvrt=n1;   ed[k].scvrt=n2; 
		ed[k].frent=i;
		ed[k].scent=-1;
		T[i][j]=k; j++;
		k++;
		}
	else
		{ed[t].scent=i;
		T[i][j]=t; j++;
		}
	t=sujm_test_fujk(ed,n1,n3,N,&s);
	if(s==-1)
		{N++;
		ed[k].frvrt=n1;   ed[k].scvrt=n3;  
		ed[k].frent=i;
		ed[k].scent=-1;
		T[i][j]=k; j++;
		k++;
		}
	else
		{ed[t].scent=i;
		T[i][j]=t; j++;
		}
	t=sujm_test_fujk(ed,n2,n3,N,&s);
	if(s==-1)
		{N++;
		ed[k].frvrt=n2;   ed[k].scvrt=n3;  
		ed[k].frent=i;
		ed[k].scent=-1;
		T[i][j]=k; j++;
		k++;
		}
	else
		{ed[t].scent=i;
		T[i][j]=t; j++;
		} 
	}
return N;
}



void tupv_fill_hagj(manif_ro *msh)
{int nel,i,**T,ned;
kt *temp;
nel=msh->e_grs;
T=(int **)malloc(nel*sizeof(int*));
for(i=0;i<nel;i++)
	T[i]=(int *)malloc(4*sizeof(int));
temp=(kt *)malloc(3*nel*sizeof(kt));
ned=vawc_dete_piln(*msh,temp,T);
msh->k_grs=ned; 
for(i=0;i<ned;i++)
	{msh->kt[i].frvrt=temp[i].frvrt;
	msh->kt[i].scvrt=temp[i].scvrt;
	msh->kt[i].frent=temp[i].frent;
	msh->kt[i].scent=temp[i].scent;
	}
for(i=0;i<nel;i++)
	{msh->entity[i].frkt=T[i][1];
	msh->entity[i].sckt=T[i][2];
	msh->entity[i].trkt=T[i][3];
	}
for(i=0;i<nel;i++)
	free(T[i]);
free(T);
free(temp);
}


void qucj_find_doct_jicw(kt e,kt *f)
{f->frent=e.frent;
f->scent=e.scent;
f->frvrt=e.frvrt;
f->scvrt=e.scvrt;
}


int wivg_node_wocp(telolf T,int w)
{int res=-1;
if(w==T.frvrt)		res=T.scvrt;
else if(w==T.scvrt)	res=T.thvrt;
else if(w==T.thvrt)	res=T.frvrt;
return res;
}


int haps_node_nepl(telolf T,int w)
{int res=-1;
if(w==T.frvrt)		res=T.thvrt;
else if(w==T.scvrt)	res=T.frvrt;
else if(w==T.thvrt)	res=T.scvrt;
return res;
}


 
void tegk_unif_zosl(manif_ro M,int e1,int e2,int pi,int pj,int pk,int *a,int *b,int *c,int *d)
{int temp;
temp=wivg_node_wocp(M.entity[e1],pi);
if(temp!=pj)
	{*a=pi; 
	*b=wivg_node_wocp(M.entity[e1],pi);
	}
else
	{*a=pj; 
	*b=wivg_node_wocp(M.entity[e1],pj);
	}
*c=wivg_node_wocp(M.entity[e1],*b);
*d=pk;
}


double homl_smal_rezt(parm A,parm B,parm C)
{double alpha,temp;
alpha=4.0;
temp=lomf_angl_serb(A,B,C);
if(temp<alpha)
	alpha=temp;
temp=lomf_angl_serb(B,C,A);
if(temp<alpha)
	alpha=temp;
temp=lomf_angl_serb(C,B,A);
if(temp<alpha)
	alpha=temp;
return alpha;
}



int pevc_find_murn_varf(parm A,parm B,parm C,parm D)
{int res,r;
double gamma,alpha;
parm *P;
gamma=1.0*MY_PI;
P=(parm *)malloc(6*sizeof(parm));
cunl_find_qedf_rewn(A,&P[0]);
cunl_find_qedf_rewn(B,&P[1]);
cunl_find_qedf_rewn(C,&P[2]);
cunl_find_qedf_rewn(D,&P[3]);
cunl_find_qedf_rewn(A,&P[4]);
cunl_find_qedf_rewn(B,&P[5]);
res=1;
for(r=1;r<=4;r++)
	{alpha=lomn_inte_cubq(P[r-1],P[r],P[r+1]);
	if(alpha>gamma)
		{res=0;
		break;
		}
	}
free(P);
return res;
}

 

int futl_take_rowg(manif_ro M,int i,int n,int m)
{int res,n1,n2,n3;
n1=M.entity[i].frvrt;
n2=M.entity[i].scvrt;
n3=M.entity[i].thvrt;
if((n==n1)&&(m==n2))res=n3;
if((n==n1)&&(m==n3))res=n2;
if((n==n2)&&(m==n3))res=n1;
if((m==n1)&&(n==n2))res=n3;
if((m==n1)&&(n==n3))res=n2;
if((m==n2)&&(n==n3))res=n1;
return res;
}



int rutl_find_senc_kamv(manif_ro M,int e)
{int res,pi,pj,pl,pk,e1,e2,a,b,c,d,tr;
double temp1,temp2,pos1,pos2;
pi=M.kt[e].frvrt;
pj=M.kt[e].scvrt;
e1=M.kt[e].frent;
e2=M.kt[e].scent;
if(e2!=-1)
	{pl=futl_take_rowg(M,e1,pi,pj);
	pk=futl_take_rowg(M,e2,pi,pj);
	tegk_unif_zosl(M,e1,e2,pi,pj,pk,&a,&b,&c,&d);
	tr=pevc_find_murn_varf(M.knot[a],M.knot[b],M.knot[c],M.knot[d]);
	if(tr==1)
		{
		temp1=homl_smal_rezt(M.knot[pi],M.knot[pj],M.knot[pk]);
		temp2=homl_smal_rezt(M.knot[pi],M.knot[pj],M.knot[pl]);
		if(temp1<temp2)
			pos1=temp1;
		else
			pos1=temp2;
		
		temp1=homl_smal_rezt(M.knot[pk],M.knot[pl],M.knot[pi]);
		temp2=homl_smal_rezt(M.knot[pk],M.knot[pl],M.knot[pj]);
		if(temp1<temp2)
			pos2=temp1;
		else
			pos2=temp2;
		
		if(pos1<pos2)
			res=0;
		else
			res=1;
		}
	else
		res=1;
	}
else
	res=1;
return res;
}



void biwg_dete_tesd(manif_ro M,int e,int m,int n,int p,int q,int *e1,int *e2,int *e3,int *e4)
{int edg1,edg2,edg3,elm1,elm2;
elm1=M.kt[e].frent;
edg1=M.entity[elm1].frkt;
edg2=M.entity[elm1].sckt;
edg3=M.entity[elm1].trkt;

if((M.kt[edg1].frvrt==m)&&(M.kt[edg1].scvrt==p))
	*e1=edg1;
if((M.kt[edg1].frvrt==p)&&(M.kt[edg1].scvrt==m))
	*e1=edg1;
if((M.kt[edg2].frvrt==m)&&(M.kt[edg2].scvrt==p))
	*e1=edg2;
if((M.kt[edg2].frvrt==p)&&(M.kt[edg2].scvrt==m))
	*e1=edg2;
if((M.kt[edg3].frvrt==m)&&(M.kt[edg3].scvrt==p))
	*e1=edg3;
if((M.kt[edg3].frvrt==p)&&(M.kt[edg3].scvrt==m))
	*e1=edg3;

if((M.kt[edg1].frvrt==n)&&(M.kt[edg1].scvrt==p))
	*e2=edg1;
if((M.kt[edg1].frvrt==p)&&(M.kt[edg1].scvrt==n))
	*e2=edg1;
if((M.kt[edg2].frvrt==n)&&(M.kt[edg2].scvrt==p))
	*e2=edg2;
if((M.kt[edg2].frvrt==p)&&(M.kt[edg2].scvrt==n))
	*e2=edg2;
if((M.kt[edg3].frvrt==n)&&(M.kt[edg3].scvrt==p))
	*e2=edg3;
if((M.kt[edg3].frvrt==p)&&(M.kt[edg3].scvrt==n))
	*e2=edg3;

elm2=M.kt[e].scent;
edg1=M.entity[elm2].frkt;
edg2=M.entity[elm2].sckt; 
edg3=M.entity[elm2].trkt;

if((M.kt[edg1].frvrt==m)&&(M.kt[edg1].scvrt==q))
	*e3=edg1;
if((M.kt[edg1].frvrt==q)&&(M.kt[edg1].scvrt==m))
	*e3=edg1;
if((M.kt[edg2].frvrt==m)&&(M.kt[edg2].scvrt==q))
	*e3=edg2;
if((M.kt[edg2].frvrt==q)&&(M.kt[edg2].scvrt==m))
	*e3=edg2;
if((M.kt[edg3].frvrt==m)&&(M.kt[edg3].scvrt==q))
	*e3=edg3;
if((M.kt[edg3].frvrt==q)&&(M.kt[edg3].scvrt==m))
	*e3=edg3;

if((M.kt[edg1].frvrt==n)&&(M.kt[edg1].scvrt==q))
	*e4=edg1;
if((M.kt[edg1].frvrt==q)&&(M.kt[edg1].scvrt==n))
	*e4=edg1;
if((M.kt[edg2].frvrt==n)&&(M.kt[edg2].scvrt==q))
	*e4=edg2;
if((M.kt[edg2].frvrt==q)&&(M.kt[edg2].scvrt==n))
	*e4=edg2;
if((M.kt[edg3].frvrt==n)&&(M.kt[edg3].scvrt==q))
	*e4=edg3;
if((M.kt[edg3].frvrt==q)&&(M.kt[edg3].scvrt==n))
	*e4=edg3;
}



void hevj_find_jerw_tefn(manif_ro *M,int e)
{int temp,m,n,p,q,elm1,elm2,e1,e2,e3,e4,sc,cs;
kt *tp;
m=M->kt[e].frvrt;
n=M->kt[e].scvrt;
elm1=M->kt[e].frent;
elm2=M->kt[e].scent;
p=futl_take_rowg(*M,elm1,n,m);
q=futl_take_rowg(*M,elm2,n,m);
sc=wivg_node_wocp(M->entity[elm1],p);
if(sc==n)	cs=1;
if(sc==m)	cs=2;
biwg_dete_tesd(*M,e,m,n,p,q,&e1,&e2,&e3,&e4);
tp=(kt *)malloc(6*sizeof(kt));
qucj_find_doct_jicw(M->kt[e1],&tp[1]);
qucj_find_doct_jicw(M->kt[e2],&tp[2]);
qucj_find_doct_jicw(M->kt[e3],&tp[3]);
qucj_find_doct_jicw(M->kt[e4],&tp[4]);
qucj_find_doct_jicw(M->kt[e],&tp[5]);

M->kt[e].frvrt=p;
M->kt[e].scvrt=q;

if(cs==1)
	{M->entity[elm1].frvrt=m;
	M->entity[elm1].scvrt=p;
	M->entity[elm1].thvrt=q;
	}
if(cs==2)
	{M->entity[elm1].frvrt=m;
	M->entity[elm1].scvrt=q;
	M->entity[elm1].thvrt=p;
	}

if(cs==1)
	{M->entity[elm2].frvrt=n;
	M->entity[elm2].scvrt=q;
	M->entity[elm2].thvrt=p;
	}
if(cs==2)
	{M->entity[elm2].frvrt=n;
	M->entity[elm2].scvrt=p;
	M->entity[elm2].thvrt=q;
	}

M->entity[elm1].frkt=e1;
M->entity[elm1].sckt=e3;
M->entity[elm1].trkt=e;

M->entity[elm2].frkt=e2;
M->entity[elm2].sckt=e4;
M->entity[elm2].trkt=e;

if(tp[1].frent!=elm1)
	temp=tp[1].frent;
else
	temp=tp[1].scent;
M->kt[e1].frent=elm1;
M->kt[e1].scent=temp;

if(tp[2].frent!=elm1)
	temp=tp[2].frent;
else
	temp=tp[2].scent;
M->kt[e2].frent=elm2;
M->kt[e2].scent=temp;

if(tp[3].frent!=elm2)
	temp=tp[3].frent;
else
	temp=tp[3].scent;
M->kt[e3].frent=elm1;
M->kt[e3].scent=temp;

if(tp[4].frent!=elm2)
	temp=tp[4].frent;
else
	temp=tp[4].scent;
M->kt[e4].frent=elm2;
M->kt[e4].scent=temp;
free(tp);
}


int cong_find_femw_jigw(manif_ro *M,int e)
{int ts,res=0;
ts=rutl_find_senc_kamv(*M,e);
if(ts==0)
	{hevj_find_jerw_tefn(M,e);
	res=1;
	}
return res;
}


int katg_find_qark_hojp(manif_ro *M)
{int res=0,ned,i,ts;
ned=M->k_grs;
for(i=0;i<ned;i++)
if(M->kt[i].scent!=-1)
	{ts=cong_find_femw_jigw(M,i);
	if(ts==1)
		res=1;
	}
return res;
}


void rodq_find_hakw_qonj(manif_ro *M,int max_leg)
{int ned,i,ts;
ned=M->k_grs;
for(i=0;i<max_leg;i++)
if(i<ned)
	{ts=katg_find_qark_hojp(M);
	if(ts==0)
		break;
	}
}

 
void jilg_conv_hekd(polygon P,mult_conn *mc)
{int i,nin;
nin=P.nb_inner_boundaries;
mc->v_grs=P.v_grs;
mc->nb_vr_outer=P.nb_local_vertices[0];
for(i=0;i<nin;i++)
	mc->nb_vr_inner[i]=P.nb_local_vertices[i+1];
for(i=0;i<P.v_grs;i++)
	{mc->vertex[i].u=P.vertex[i].u;
	mc->vertex[i].v=P.vertex[i].v;
	}
mc->nb_inner_polygons=nin;
}



int delaunay_triangulate(polygon P,manif_ro *msh,int *forc_term)
{int nin,n,suc,max_leg=20,f_trm;
mult_conn mc;
*forc_term=0;
nin=P.nb_inner_boundaries;
n=P.v_grs;
mc.vertex=(parm *)malloc(n*sizeof(parm));
mc.zt=(double *)malloc(n*sizeof(double));
mc.flag=(int *)malloc(n*sizeof(int));
mc.nb_vr_inner=(int *)malloc(nin*sizeof(int));
mc.mapglob=(int *)malloc(n*sizeof(int));
jilg_conv_hekd(P,&mc);
suc=vorg_find_qach_nujt(mc,msh,&f_trm);
if(f_trm==1)
	{*forc_term=1;
	fprintf(tmpout,"force term: vorg_find_qach_nujt() in delaunay_triangulate()\n");
	free(mc.nb_vr_inner);
	free(mc.vertex);
	free(mc.zt);
	free(mc.flag);
	free(mc.mapglob);
	return FAILURE;
	}
if(suc==FAILURE)
	{msh->e_grs=0;
	msh->n_grs=0;
	}
else
	{tupv_fill_hagj(msh);
	rodq_find_hakw_qonj(msh,max_leg);
	}
free(mc.nb_vr_inner);
free(mc.vertex);
free(mc.zt);
free(mc.flag);
free(mc.mapglob);
return suc;
} 


