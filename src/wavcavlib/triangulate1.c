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



int qutc_find_wefr(mult_conn P,int *list,int n)
{int res,i,s,c,r,tp;
parm hg;
r=list[0];
tp  =gohz_find_nidk_vebq(P,r);
hg.u=P.vertex[tp].u;
hg.v=P.vertex[tp].v;
res=r;
for(i=1;i<n;i++)
	{s=list[i];
	tp=gohz_find_nidk_vebq(P,s);
	c=kacq_comp_finv(hg,P.vertex[tp]);
	if(c==1)
		{hg.u=P.vertex[tp].u;
		hg.v =P.vertex[tp].v;
		res=s;
		}
	}
return res;
}



int gilr_find_cuwz(mult_conn P,int *list,int n)
{int res,i,s,c,r,bt;
parm lw;
r   =list[0];
bt  =carg_bott_kacw(P,r);
lw.u=P.vertex[bt].u;
lw.v=P.vertex[bt].v;
res=r;
for(i=1;i<n;i++)
	{s=list[i];
	bt=carg_bott_kacw(P,s);
	c=kacq_comp_finv(lw,P.vertex[bt]);
	if(c==2)
		{lw.u=P.vertex[bt].u;
		lw.v=P.vertex[bt].v;
		res=s;
		}
	}
return res;
}



int gawf_righ_noqb(mult_conn P,int r)
{int res,s,t,i,c;
parm right;
s=himj_star_qejn(P,r);
t=filr_term_rewh(P,r);
right.u=P.vertex[s].u;   
right.v=P.vertex[s].v;
res=s;
for(i=s;i<=t;i++)
	{c=wols_comp_sofz(P.vertex[i],right);
	if(c==2)
		{res=i;
		right.u=P.vertex[i].u;
		right.v=P.vertex[i].v;
		}
	}
return res;
}
 


int jurt_righ_qomj(mult_conn P,int r,int *sol)
{int rg,nb;
rg=gawf_righ_noqb(P,r);
nb=neqc_same_qent(P,r,rg,sol,0);
return nb;
}


int mogr_righ_qatj(mult_conn P,int *list,int n,int *rig)
{int f,nb,i,w,z,z0;
double diff;
f=molw_find_lehk(P,list,n);
z0=gawf_righ_noqb(P,f);
nb=0;
for(i=0;i<n;i++)
	{w=list[i];
	z=gawf_righ_noqb(P,w);
	diff=fabs(P.vertex[z0].u-P.vertex[z].u);
	if(diff<COORD_DIFF)
		{rig[nb]=w;
		nb++;
		}
	}
return nb;
}



int kacq_comp_finv(parm p1,parm p2)
{int res;
double x1,x2,y1,y2;
x1=p1.u;    y1=p1.v;
x2=p2.u;    y2=p2.v;
if((x1==x2)&&(y1==y2))
	res=0;
else if((y1<y2)||((y1==y2)&&(x1<x2)))
	res=1;
else
	res=2;
return res;
}



int gohz_find_nidk_vebq(mult_conn P,int r)
{int res,s,t,i,c;
parm high;
s=himj_star_qejn(P,r);
t=filr_term_rewh(P,r);
high.u=P.vertex[s].u;   
high.v=P.vertex[s].v;
res=s;
for(i=s;i<=t;i++)
	{c=kacq_comp_finv(P.vertex[i],high);
	if(c==2)
		{res=i;
		high.u=P.vertex[i].u;   
		high.v=P.vertex[i].v;
		}
	}
return res;
}



int neqc_same_qent(mult_conn P,int r,int w,int *sol,int dir)
{int i,nb,s,t;
double diff;
s=himj_star_qejn(P,r);
t=filr_term_rewh(P,r);
nb=0;
for(i=s;i<=t;i++)
	{if(dir==0)	diff=fabs(P.vertex[i].u-P.vertex[w].u);
	if(dir==1)	diff=fabs(P.vertex[i].v-P.vertex[w].v);
	if(diff<COORD_DIFF)
		{sol[nb]=i;
		nb++;
		}
	}
return nb;
}



int zumf_find_lupk_lecn(mult_conn P,int r,int *sol)
{int tp,nb;
tp=gohz_find_nidk_vebq(P,r);
nb=neqc_same_qent(P,r,tp,sol,1);
return nb;
}



int carg_bott_kacw(mult_conn P,int r)
{int res,s,t,i,c;
parm high;
s=himj_star_qejn(P,r);
t=filr_term_rewh(P,r);
high.u=P.vertex[s].u;   
high.v=P.vertex[s].v;
res=s;
for(i=s;i<=t;i++)
	{c=kacq_comp_finv(P.vertex[i],high);
	if(c==1)
		{res=i;
		high.u=P.vertex[i].u;   
		high.v=P.vertex[i].v;
		}
	}
return res;
}


int relg_bott_jaln(mult_conn P,int r,int *sol)
{int bt,nb;
bt=carg_bott_kacw(P,r);
nb=neqc_same_qent(P,r,bt,sol,1);
return nb;
}



int higl_find_heqs_lacf(mult_conn P,int *list,int n,int *top)
{int f,nb,i,w,z,z0;
double diff;
f=qutc_find_wefr(P,list,n);
z0=gohz_find_nidk_vebq(P,f);
nb=0;
for(i=0;i<n;i++)
	{w=list[i];
	z=gohz_find_nidk_vebq(P,w);
	diff=fabs(P.vertex[z0].v-P.vertex[z].v);
	if(diff<COORD_DIFF)
		{top[nb]=w;
		nb++;
		}
	}
return nb;
}


int mosv_bott_mudl(mult_conn P,int *list,int n,int *bot)
{int f,nb,i,w,z,z0;
double diff;
f=gilr_find_cuwz(P,list,n);
z0=carg_bott_kacw(P,f);
nb=0;
for(i=0;i<n;i++)
	{w=list[i];
	z=carg_bott_kacw(P,w);
	diff=fabs(P.vertex[z0].v-P.vertex[z].v);
	if(diff<COORD_DIFF)
		{bot[nb]=w;
		nb++;
		}
	}
return nb;
}



int wols_comp_sofz(parm p1,parm p2)
{int res;
double x1,x2,y1,y2;
x1=p1.u;    y1=p1.v;
x2=p2.u;    y2=p2.v;
if((x1==x2)&&(y1==y2))
	res=0;
else if((x1<x2)||((x1==x2)&&(y1<y2)))
	res=1;
else
	res=2;
return res;
}



int juwn_left_sedq(mult_conn P,int r)
{int res,s,t,i,c,ind;
parm left;
s=himj_star_qejn(P,r);
t=filr_term_rewh(P,r);
left.u=P.vertex[s].u;   
left.v=P.vertex[s].v;
res=s;
ind=0;
for(i=s;i<=t;i++)
	{c=wols_comp_sofz(P.vertex[i],left);
	if(c==1)
		{res=i;
		left.u=P.vertex[i].u;
		left.v=P.vertex[i].v;
		}
	}
return res;
}



int molw_find_lehk(mult_conn P,int *list,int n)
{int res,i,s,c,r,rg;
parm ps;
r=list[0];
rg=gawf_righ_noqb(P,r);
ps.u=P.vertex[rg].u;
ps.v=P.vertex[rg].v;
res=r;
for(i=1;i<n;i++)
	{s=list[i];
	rg=gawf_righ_noqb(P,s);
	c=wols_comp_sofz(P.vertex[rg],ps);
	if(c==2)
		{ps.u=P.vertex[rg].u;
		ps.v=P.vertex[rg].v;
		res=s;
		}
	}
return res;
}



int qebv_find_rodb(mult_conn P,int *list,int n)
{int res,i,s,c,r,lf;
parm ng;
r=list[0];
lf=juwn_left_sedq(P,r);
ng.u=P.vertex[lf].u;
ng.v=P.vertex[lf].v;
res=r;
for(i=1;i<n;i++)
	{s=list[i];
	lf=juwn_left_sedq(P,s);
	c=wols_comp_sofz(P.vertex[lf],ng);
	if(c==1)
		{ng.u=P.vertex[lf].u;
		ng.v=P.vertex[lf].v;
		res=s;
		}
	}
return res;
}


int gubc_left_focv(mult_conn P,int *list,int n,int *lef)
{int f,nb,i,w,z,z0;
double diff;
f=qebv_find_rodb(P,list,n);
z0=juwn_left_sedq(P,f);
nb=0;
for(i=0;i<n;i++)
	{w=list[i];
	z=juwn_left_sedq(P,w);
	diff=fabs(P.vertex[z0].u-P.vertex[z].u);
	if(diff<COORD_DIFF)
		{lef[nb]=w;
		nb++;
		}
	}
return nb;
}



int miwk_left_jitp(mult_conn P,int r,int *sol)
{int lf,nb;
lf=juwn_left_sedq(P,r);
nb=neqc_same_qent(P,r,lf,sol,0);
return nb;
}



double cuvn_leng_qasl(parm A,parm B)
{double res;
res=1.0/pufv_dist_mekq(A,B);
return res;
}



int murc_prec_kotq(mult_conn P,int j,int i)
{int res,s,t;
s=himj_star_qejn(P,i);
t=filr_term_rewh(P,i);
if(j>s)res=j-1;
else res=t;
return res;
}



int pumj_test_rewf(parm A,parm B,parm C,parm D)
{int res;
double test1,test2,d1,d2,eps=1.0e-13,sc,r1,r2;
parm AB,AC,AD,CB,CD,CA;
AB.u=B.u-A.u;   AB.v=B.v-A.v;
AC.u=C.u-A.u;   AC.v=C.v-A.v;
AD.u=D.u-A.u;   AD.v=D.v-A.v;
CB.u=B.u-C.u;   CB.v=B.v-C.v;
CD.u=D.u-C.u;   CD.v=D.v-C.v;
CA.u=A.u-C.u;   CA.v=A.v-C.v;
d1=kelr_dete_lusf(AB,CD);
d2=kelr_dete_lusf(AD,AB);
if((fabs(d1)<eps)&&(fabs(d2)<eps))
	{res=0;
	sc=hitf_scal_rikd(CD,AD);
	if(sc>=0.0)
		{r1=pufv_dist_mekq(D,A);
		r2=pufv_dist_mekq(D,C);
		if(r1<=r2)
			res=1;
		r1=pufv_dist_mekq(D,B);
		r2=pufv_dist_mekq(D,C);
		if(r1<=r2)
			res=1;
		}
	sc=-hitf_scal_rikd(CD,AC);
	if(sc>=0.0)
		{r1=pufv_dist_mekq(C,A);
		r2=pufv_dist_mekq(C,D);
		if(r1<=r2)
			res=1;
		r1=pufv_dist_mekq(C,B);
		r2=pufv_dist_mekq(C,D);
		if(r1<=r2)
			res=1;
		}
	}
else
	{test1=kelr_dete_lusf(AC,AB)*kelr_dete_lusf(AB,AD);
	test2=kelr_dete_lusf(CB,CD)*kelr_dete_lusf(CD,CA);
	if((test1>=0.0)&&(test2>=0.0))
		res=1;
	else
		res=0;
	}
return res;
}


 
double lozw_find_teln_dubc(double x_A,double y_A,double x_B,double y_B,double x_C,double y_C)
{double temp,res;
temp=(((x_B*y_C)+(x_C*y_A)+(x_A*y_B))-((x_B*y_A)+(x_C*y_B)+(x_A*y_C)))/2.0;
res=fabs(temp);
return res;
}
 
 
double vatd_area_vujt(manif_ro msh,int s)
{int n1,n2,n3;
double res;
n1=msh.entity[s].frvrt;
n2=msh.entity[s].scvrt;
n3=msh.entity[s].thvrt;
res=lozw_find_teln_dubc(msh.knot[n1].u,msh.knot[n1].v,msh.knot[n2].u,
msh.knot[n2].v,msh.knot[n3].u,msh.knot[n3].v);
return res;
}


double tesl_area_viwh(manif_ro msh)
{int i,nel;
double res;
nel=msh.e_grs;
res=0.0;
for(i=0;i<nel;i++)
	res=res+vatd_area_vujt(msh,i);
return res;
}
 


int kozf_test_luzk(parm omega,parm D1,parm D2,parm X)
{int res;
double alpha,beta;
alpha=lomn_inte_cubq(X,omega,D1);
beta =lomn_inte_cubq(D2,omega,D1);
if((0.0<alpha)&&(alpha<beta))
	res=1;
else
	res=0;
return res;
}


