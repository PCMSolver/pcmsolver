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
#include <stdlib.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "triang.h"


void jonc_find_qifn_fupw(manif_ro in,manif_ro *out)
{int i,nnd,nel,ned;
nnd=in.n_grs;
for(i=0;i<nnd;i++)
	{out->knot[i].u=in.knot[i].u;
	out->knot[i].v=in.knot[i].v;
	}
out->n_grs=nnd;

nel=in.e_grs;
for(i=0;i<nel;i++)
	{out->entity[i].frvrt=in.entity[i].frvrt;
	out->entity[i].scvrt=in.entity[i].scvrt;
	out->entity[i].thvrt=in.entity[i].thvrt;
	out->entity[i].frkt=in.entity[i].frkt;
	out->entity[i].sckt=in.entity[i].sckt;
	out->entity[i].trkt=in.entity[i].trkt;
	}
out->e_grs=nel;

ned=in.k_grs;
for(i=0;i<ned;i++)
	{out->kt[i].frent=in.kt[i].frent;
	out->kt[i].scent=in.kt[i].scent;
	out->kt[i].frvrt=in.kt[i].frvrt;
	out->kt[i].scvrt=in.kt[i].scvrt;
	}
out->k_grs=ned;
}



double purt_aspe_kudl(manif_ro msh,int s)
{int i,j,k;
double res,*r,temp;
parm A,B,C;
i=msh.entity[s].frvrt;
j=msh.entity[s].scvrt;
k=msh.entity[s].thvrt;
A.u=msh.knot[i].u;   A.v=msh.knot[i].v;
B.u=msh.knot[j].u;   B.v=msh.knot[j].v;
C.u=msh.knot[k].u;   C.v=msh.knot[k].v;
r=(double *)malloc(3*sizeof(double));
r[0]=lomf_angl_serb(A,B,C);
r[1]=lomf_angl_serb(B,A,C);
r[2]=lomf_angl_serb(C,A,B);
temp=4.0;
for(i=0;i<3;i++)
	if(r[i]<temp)
		temp=r[i];
free(r);
res=1.0/temp;
return res;
}



int vokj_long_vern(manif_ro msh,int s)
{int i,j,k,res,*e;
double *r,temp;
e=(int *)malloc(3*sizeof(int));
e[0]=msh.entity[s].frkt;
e[1]=msh.entity[s].sckt;
e[2]=msh.entity[s].trkt;
r=(double *)malloc(3*sizeof(double));
for(i=0;i<3;i++)
	{j=msh.kt[e[i]].frvrt;
	k=msh.kt[e[i]].scvrt;
	r[i]=pufv_dist_mekq(msh.knot[j],msh.knot[k]);
	}
temp=0.0;
for(i=0;i<3;i++)
	if(temp<r[i])
		{temp=r[i];
		res=e[i];
		}
free(r);
free(e);
return res;
}



int jokq_test_julg(parm A,parm B,parm C,parm P)
{int res,ts;
double lm1,lm2,lm3;
ts=cerv_segm_gusk(A,B,P);
if(ts==1)	return 1;
ts=cerv_segm_gusk(B,C,P);
if(ts==1)	return 1;
ts=cerv_segm_gusk(A,C,P);
if(ts==1)	return 1;
lm1=vupq_lamb_qofc(A.u,A.v,B.u,B.v,C.u,C.v,P.u,P.v);
lm2=dopg_lamb_nupd(A.u,A.v,B.u,B.v,C.u,C.v,P.u,P.v);
lm3=mofr_lamb_powg(A.u,A.v,B.u,B.v,C.u,C.v,P.u,P.v);
if((lm1>=0.0)&&(lm2>=0.0)&&(lm3>=0.0))
	res=1;
else	
	res=0;
return res;
}



int fegb_dete_bejw(manif_ro M,int e)
{int res,pi,pj,pl,pk,e1,e2,cq,cr,cs,ct;
double temp1,temp2,pos1,pos2;
pi=M.kt[e].frvrt;
pj=M.kt[e].scvrt;
e1=M.kt[e].frent;
e2=M.kt[e].scent;
if(e2!=-1)
	{pl=futl_take_rowg(M,e1,pi,pj);
	pk=futl_take_rowg(M,e2,pi,pj);
	cq=jokq_test_julg(M.knot[pi],M.knot[pk],M.knot[pl],M.knot[pj]);
	cr=jokq_test_julg(M.knot[pi],M.knot[pj],M.knot[pl],M.knot[pk]);
	cs=jokq_test_julg(M.knot[pi],M.knot[pk],M.knot[pj],M.knot[pl]);
	ct=jokq_test_julg(M.knot[pj],M.knot[pk],M.knot[pl],M.knot[pi]);
	if((cq==0)&&(cr==0)&&(cs==0)&&(ct==0))
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



int lers_find_judv(manif_ro M,int i,int r)
{int res,e1,e2,e3;
e1=M.entity[i].frkt;
e2=M.entity[i].sckt;
e3=M.entity[i].trkt;
if((M.kt[e1].frvrt!=r)&&(M.kt[e1].scvrt!=r))
	res=e1;
if((M.kt[e2].frvrt!=r)&&(M.kt[e2].scvrt!=r))
	res=e2;
if((M.kt[e3].frvrt!=r)&&(M.kt[e3].scvrt!=r))
	res=e3;
return res;
} 


void lodv_lega_holt(manif_ro *M,int e,int r)
{int temp,elm1,elm2,f;
temp=fegb_dete_bejw(*M,e);
if(temp==0)
	{hevj_find_jerw_tefn(M,e);
	elm1=M->kt[e].frent;
	elm2=M->kt[e].scent;
	f=lers_find_judv(*M,elm1,r);
	lodv_lega_holt(M,f,r);
	f=lers_find_judv(*M,elm2,r);
	lodv_lega_holt(M,f,r);
	}
}



int rehs_dete_pevt(manif_ro M,int s,int ni,int nj)
{int res,e1,e2,e3;
e1=M.entity[s].frkt;
e2=M.entity[s].sckt;
e3=M.entity[s].trkt;
if((M.kt[e1].frvrt==ni)&&(M.kt[e1].scvrt==nj))
	res=e1;
if((M.kt[e1].frvrt==nj)&&(M.kt[e1].scvrt==ni))
	res=e1;
if((M.kt[e2].frvrt==ni)&&(M.kt[e2].scvrt==nj))
	res=e2;
if((M.kt[e2].frvrt==nj)&&(M.kt[e2].scvrt==ni))
	res=e2;
if((M.kt[e3].frvrt==ni)&&(M.kt[e3].scvrt==nj))
	res=e3;
if((M.kt[e3].frvrt==nj)&&(M.kt[e3].scvrt==ni))
	res=e3;
return res;
}



void fikg_inse_cenp(manif_ro *M,parm P,int s,int ni,int nj)
{int ned,nel,nnd,nk,nl,t,e,e1,e2,e3,e4;

nk=futl_take_rowg(*M,s,ni,nj);
e=rehs_dete_pevt(*M,s,ni,nj);
if(M->kt[e].frent==s)
	t=M->kt[e].scent;
else
	t=M->kt[e].frent;
nl=futl_take_rowg(*M,t,ni,nj);

e1=rehs_dete_pevt(*M,s,nk,nj);
e2=rehs_dete_pevt(*M,t,nj,nl);
e3=rehs_dete_pevt(*M,s,nk,ni);
e4=rehs_dete_pevt(*M,t,ni,nl);
ned=M->k_grs;
nel=M->e_grs;
nnd=M->n_grs;

M->knot[nnd].u=P.u;
M->knot[nnd].v=P.v;
M->n_grs=nnd+1;

M->kt[ned].frvrt=nnd;    M->kt[ned].scvrt=nk;
M->kt[ned].frent=s;      M->kt[ned].scent=nel;
M->kt[ned+1].frvrt=nnd;  M->kt[ned+1].scvrt=ni;
M->kt[ned+1].frent=nel;  M->kt[ned+1].scent=nel+1;
M->kt[ned+2].frvrt=nnd;  M->kt[ned+2].scvrt=nl;
M->kt[ned+2].frent=t;    M->kt[ned+2].scent=nel+1;
M->k_grs=ned+3;

M->entity[nel].frvrt=ni;
M->entity[nel].scvrt=nnd;
M->entity[nel].thvrt=nk;
M->entity[nel].frkt=e3;
M->entity[nel].sckt=ned;
M->entity[nel].trkt=ned+1;
M->entity[nel+1].frvrt=ni;
M->entity[nel+1].scvrt=nnd;
M->entity[nel+1].thvrt=nl;
M->entity[nel+1].frkt=e4;
M->entity[nel+1].sckt=ned+1;
M->entity[nel+1].trkt=ned+2;
M->e_grs=nel+2;

M->kt[e].frvrt=nj;
M->kt[e].scvrt=nnd;

M->entity[s].frvrt=nj;
M->entity[s].scvrt=nnd;
M->entity[s].thvrt=nk;

M->entity[t].frvrt=nj;
M->entity[t].scvrt=nnd;
M->entity[t].thvrt=nl;

M->entity[s].frkt=e1;
M->entity[s].sckt=e;
M->entity[s].trkt=ned;

M->entity[t].frkt=e2;
M->entity[t].sckt=e;
M->entity[t].trkt=ned+2;

if(M->kt[e3].frent==s)M->kt[e3].frent=nel;
if(M->kt[e3].scent==s)M->kt[e3].scent=nel;

if(M->kt[e4].frent==t)M->kt[e4].frent=nel+1;
if(M->kt[e4].scent==t)M->kt[e4].scent=nel+1;
lodv_lega_holt(M,e1,nnd);
lodv_lega_holt(M,e2,nnd);
lodv_lega_holt(M,e3,nnd);
lodv_lega_holt(M,e4,nnd);
}



void padl_inse_qefr(manif_ro *M,parm P,int s,int ni,int nj)
{int ned,nel,nnd,nk,e,e1,e2;

nk=futl_take_rowg(*M,s,ni,nj);
e=rehs_dete_pevt(*M,s,ni,nj);

e1=rehs_dete_pevt(*M,s,nk,nj);
e2=rehs_dete_pevt(*M,s,nk,ni);
ned=M->k_grs;
nel=M->e_grs;
nnd=M->n_grs;

M->knot[nnd].u=P.u;
M->knot[nnd].v=P.v;
M->n_grs=nnd+1;

M->kt[ned].frvrt=nnd;    M->kt[ned].scvrt=nk;
M->kt[ned].frent=s;      M->kt[ned].scent=nel;
M->kt[ned+1].frvrt=nnd;  M->kt[ned+1].scvrt=ni;
M->kt[ned+1].frent=nel;  M->kt[ned+1].scent=-1;
M->k_grs=ned+2;

M->entity[nel].frvrt=ni;
M->entity[nel].scvrt=nnd;
M->entity[nel].thvrt=nk;
M->entity[nel].frkt=e2;
M->entity[nel].sckt=ned;
M->entity[nel].trkt=ned+1;
M->e_grs=nel+1;

M->kt[e].frvrt=nj;
M->kt[e].scvrt=nnd;

M->entity[s].frvrt=nj;
M->entity[s].scvrt=nnd;
M->entity[s].thvrt=nk;

M->entity[s].frkt=e1;
M->entity[s].sckt=e;
M->entity[s].trkt=ned;

if(M->kt[e2].frent==s)M->kt[e2].frent=nel;
if(M->kt[e2].scent==s)M->kt[e2].scent=nel;
lodv_lega_holt(M,e1,nnd);
lodv_lega_holt(M,e2,nnd);
}


void sojv_inse_rost(manif_ro *M,parm P,int s,int ni,int nj)
{int e;
e=rehs_dete_pevt(*M,s,ni,nj);
if(M->kt[e].scent==-1)
	padl_inse_qefr(M,P,s,ni,nj);
else
	fikg_inse_cenp(M,P,s,ni,nj);
}


int qukr_comm_zoms(manif_ro M,int i,int j)
{int res,a,b,c,d;
a=M.kt[i].frvrt;
b=M.kt[i].scvrt;
c=M.kt[j].frvrt;
d=M.kt[j].scvrt;
if(a==c)   res=a;
if(b==c)   res=b;
if(a==d)   res=a;
if(b==d)   res=b;
return res;
}



void lask_inse_gifk(manif_ro *M,parm P,int s)
{int ned,nel,nnd,ni,nj,nk,e1,e2,e3;
e1=M->entity[s].frkt;
e2=M->entity[s].sckt;
e3=M->entity[s].trkt;
ni=qukr_comm_zoms(*M,e1,e2);
nj=qukr_comm_zoms(*M,e1,e3);
nk=qukr_comm_zoms(*M,e2,e3);
ned=M->k_grs;
nel=M->e_grs;
nnd=M->n_grs;

M->knot[nnd].u=P.u;
M->knot[nnd].v=P.v;
M->n_grs=nnd+1;

M->kt[ned].frvrt=ni;      M->kt[ned].scvrt=nnd;
M->kt[ned].frent=s;       M->kt[ned].scent=nel;
M->kt[ned+1].frvrt=nj;    M->kt[ned+1].scvrt=nnd;
M->kt[ned+1].frent=s;     M->kt[ned+1].scent=nel+1;
M->kt[ned+2].frvrt=nk;    M->kt[ned+2].scvrt=nnd;
M->kt[ned+2].frent=nel;   M->kt[ned+2].scent=nel+1;
M->k_grs=ned+3;

M->entity[nel].frvrt=nnd;
M->entity[nel].scvrt=ni;
M->entity[nel].thvrt=nk;
M->entity[nel].frkt=e2;
M->entity[nel].sckt=ned;
M->entity[nel].trkt=ned+2;
M->entity[nel+1].frvrt=nnd;
M->entity[nel+1].scvrt=nj;
M->entity[nel+1].thvrt=nk;
M->entity[nel+1].frkt=e3;
M->entity[nel+1].sckt=ned+1;
M->entity[nel+1].trkt=ned+2;
M->e_grs=nel+2;

M->entity[s].frvrt=nnd;
M->entity[s].scvrt=ni;
M->entity[s].thvrt=nj;
M->entity[s].frkt=e1;
M->entity[s].sckt=ned;
M->entity[s].trkt=ned+1;

;

if(M->kt[e2].frent==s)M->kt[e2].frent=nel;
if(M->kt[e2].scent==s)M->kt[e2].scent=nel;

if(M->kt[e3].frent==s)M->kt[e3].frent=nel+1;
if(M->kt[e3].scent==s)M->kt[e3].scent=nel+1;
lodv_lega_holt(M,e1,nnd);
lodv_lega_holt(M,e2,nnd);
lodv_lega_holt(M,e3,nnd);
}



void fapn_sele_holp(manif_ro mshin,double anisotropy,double accuracy,manif_ro *mshout)
{int i,nel,ned,nnd,n1,n2,n3,lged,adj,lged2,s;
double er,ar,ar2;
manif_ro msh;
parm P;
nel=mshin.e_grs;
nnd=mshin.n_grs;
ned=mshin.k_grs;
msh.entity=(telolf *)malloc(3*nel*sizeof(telolf));
msh.knot=(parm *)malloc((nnd+nel)*sizeof(parm));
msh.kt=(kt *)malloc((nnd+4*nel+20)*sizeof(kt));
jonc_find_qifn_fupw(mshin,&msh);
for(i=0;i<nel;i++)
	{ar=purt_aspe_kudl(msh,i);
	if(ar>=anisotropy)
		{lged=vokj_long_vern(msh,i);
		if(msh.kt[lged].frent==i)
			adj=msh.kt[lged].scent;
		else
			adj=msh.kt[lged].frent;
		ar2=0.0;
		if(adj!=-1)
			{ar2=purt_aspe_kudl(msh,adj);
			lged2=vokj_long_vern(msh,adj);
			}
		if(((ar2>=anisotropy)&&(lged==lged2))||(adj==-1))
			{n1=msh.kt[lged].frvrt;
			n2=msh.kt[lged].scvrt;
			if(msh.kt[lged].scent!=-1)
				{P.u=0.5*msh.knot[n1].u+0.5*msh.knot[n2].u;
				P.v=0.5*msh.knot[n1].v+0.5*msh.knot[n2].v;
				s=msh.kt[lged].frent;
				sojv_inse_rost(&msh,P,s,n1,n2);
				}
			}
		}
	else
		{er=vatd_area_vujt(msh,i);
		if(er>=accuracy)
			{n1=msh.entity[i].frvrt;
			n2=msh.entity[i].scvrt;
			n3=msh.entity[i].thvrt;
			P.u=(msh.knot[n1].u+msh.knot[n2].u+msh.knot[n3].u)/3.0;
			P.v=(msh.knot[n1].v+msh.knot[n2].v+msh.knot[n3].v)/3.0;
			s=i;
			lask_inse_gifk(&msh,P,s);
			}
		}
	}
jonc_find_qifn_fupw(msh,mshout);
free(msh.entity);
free(msh.kt);
free(msh.knot);
}



void multiple_refinement(manif_ro mshin,int max,double anisotropy,
double accuracy,manif_ro *mshout)
{int i,nel,nnd,ned;
manif_ro tempin,tempout;
nel=mshin.e_grs;
nnd=mshin.n_grs;
ned=mshin.k_grs;
tempin.entity=(telolf *)malloc(nel*sizeof(telolf));
tempin.knot=(parm *)malloc(nnd*sizeof(parm));
tempin.kt=(kt *)malloc(ned*sizeof(kt));
jonc_find_qifn_fupw(mshin,&tempin);
for(i=0;i<max;i++)
	{nel=tempin.e_grs;
	nnd=tempin.n_grs;
	ned=tempin.k_grs;
	tempout.entity=(telolf *)malloc(3*nel*sizeof(telolf));
	tempout.knot=(parm *)malloc((nnd+nel)*sizeof(parm));
	tempout.kt=(kt *)malloc((nnd+4*nel+20)*sizeof(kt));
	fapn_sele_holp(tempin,anisotropy,accuracy,&tempout);
	
	free(tempin.kt);
	free(tempin.entity);
	free(tempin.knot);
	
	if(i<max-1)
		{nel=tempout.e_grs;
		nnd=tempout.n_grs;
		ned=tempout.k_grs;
		tempin.entity=(telolf *)malloc(nel*sizeof(telolf));
		tempin.knot=(parm *)malloc(nnd*sizeof(parm));
		tempin.kt=(kt *)malloc(ned*sizeof(kt));
		jonc_find_qifn_fupw(tempout,&tempin);
		}
	else
		jonc_find_qifn_fupw(tempout,mshout);
	
	free(tempout.kt);
	free(tempout.entity);
	free(tempout.knot);
	}
}


