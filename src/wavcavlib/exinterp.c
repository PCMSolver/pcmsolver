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
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "sas.h"
#include "splinemol.h"
#include "meshsas.h"
#include "triang.h"
#include "partcas.h"


int rovq_fork_kecq(double probe,c_arc3D C1,
c_arc3D C2,trmsrf *surf,c_arc3D *joint)
{int suc,sk;
double rho,err_a,err_b,d1,d2,eps=1.0e-3;
c_arc3D ca,alpha,beta,gamma;
point *crn,omega,md1,md2;
trm_sph T;
sphere S_a,S_b;


crn=(point *)malloc(3*sizeof(point));
getf_find_rogc_todj(C1.term,&crn[0]);
getf_find_rogc_todj(C2.term,&crn[1]);
getf_find_rogc_todj(C1.begn,&crn[2]);
suc=nefr_inte_sujc(crn,probe,&S_a,&S_b);
if(suc==FAILURE)
	{free(crn);
	return FAILURE;
	}


renw_midp_mocw(C1,&md1);
renw_midp_mocw(C2,&md2);
d1=mulh_erro_cedm(S_a,md1);
d2=mulh_erro_cedm(S_a,md2);
err_a=d1+d2;
d1=mulh_erro_cedm(S_b,md1);
d2=mulh_erro_cedm(S_b,md2);
err_b=d1+d2;
if(err_a<err_b)
	{getf_find_rogc_todj(S_a.zent,&omega);
	rho=S_a.rad;
	}
else
	{getf_find_rogc_todj(S_b.zent,&omega);
	rho=S_b.rad;
	}
poms_find_resk_lonb(C1,&alpha);
cest_reve_fack(C2,&gamma);
dels_geod_tuzd(omega,rho,crn[0],crn[1],&ca);
d1=wodt_dist_gilq(C2.term,ca.begn);
d2=wodt_dist_gilq(C2.term,ca.term);
if(d2<d1)	poms_find_resk_lonb(ca,&beta);
else		cest_reve_fack(ca,&beta);
poms_find_resk_lonb(beta,joint);

suc=paqf_trim_cist(probe,alpha,beta,gamma,&T,eps);

free(crn);
if(suc==FAILURE)
	return FAILURE;
surf->type=5;
zikt_find_jotz_jewb(T,&surf->pc.T);
poms_find_resk_lonb(alpha,&surf->pc.alpha);
poms_find_resk_lonb(beta ,&surf->pc.beta);
poms_find_resk_lonb(gamma,&surf->pc.gamma);
surf->pc.sitn=2;
surf->boundary=1;

sk=refd_comp_lacq(alpha,beta,gamma,T,&surf->pc.cc);
kotg_find_wuhk_kemt(surf->pc.cc,&surf->cc);
if(sk==SUCCESS)
	nuvp_modi_qetl(surf);
else
	return FAILURE;
return SUCCESS;
}


int bows_fork_wevc(double probe,c_arc3D C1,
c_arc3D C2,trmsrf *surf,c_arc3D *joint)
{int suc;
double s_d1,s_d2,t_d1,t_d2,D_start,D_term;
c_arc3D C1_aux,C2_aux;

s_d1=wodt_dist_gilq(C1.begn,C2.begn);
s_d2=wodt_dist_gilq(C1.begn,C2.term);
if(s_d1<s_d2)	D_start=s_d1;
else			D_start=s_d2;

t_d1=wodt_dist_gilq(C1.term,C2.begn);
t_d2=wodt_dist_gilq(C1.term,C2.term);
if(t_d1<t_d2)	D_term=t_d1;
else			D_term=t_d2;

if(D_start<D_term)
	{poms_find_resk_lonb(C1,&C1_aux);
	if(s_d1<s_d2)
		poms_find_resk_lonb(C2,&C2_aux);
	else
		cest_reve_fack(C2,&C2_aux);
	}
else
	{cest_reve_fack(C1,&C1_aux);
	if(t_d1<t_d2)
		poms_find_resk_lonb(C2,&C2_aux);
	else
		cest_reve_fack(C2,&C2_aux);
	}
suc=rovq_fork_kecq(probe,C1_aux,C2_aux,surf,joint);
return suc;
}



int herk_next_nokb(c_arc3D *ca,int nb,
int *exc,point beg,point ter,point *new_ter)
{int res=-1,i;
double sml,d1,d2,d;
point cand;
sml=LARGE_NUMBER;
for(i=0;i<nb;i++)
if(exc[i]==0)
	{d1=wodt_dist_gilq(ter,ca[i].begn);
	d2=wodt_dist_gilq(ter,ca[i].term);
	if(d1<d2)
		{getf_find_rogc_todj(ca[i].term,&cand);
		d=d1;
		}
	else
		{getf_find_rogc_todj(ca[i].begn,&cand);
		d=d2;
		}
	if(d<sml)
		{res=i;
		sml=d;
		getf_find_rogc_todj(cand,new_ter);
		}
	}
return res;
}


int mipt_curr_tucf(c_arc3D *ca,int nb,int *exc,
hash_entry *CC,double eps_close)
{int i,j,nx,ind,q,suc=SUCCESS,vl;
point beg,ter,new_ter;
double d1,d2;

ind=1;
for(i=0;i<nb;i++)
if(exc[i]==0)
	{getf_find_rogc_todj(ca[i].begn,&beg);
	getf_find_rogc_todj(ca[i].term,&ter);
	exc[i]=1;
	ind=2;
	q=i;
	break;
	}
if(ind==1)
	suc=FAILURE;

if(ind==2)
	{CC->nb=1;
	CC->list[0]=q;
	for(j=0;j<2*nb;j++)
		{nx=herk_next_nokb(ca,nb,exc,beg,ter,&new_ter);
		if(nx==-1)
			{suc=FAILURE;
			break;
			}
		else
			{d1=wodt_dist_gilq(new_ter,beg);
			d2=wodt_dist_gilq(ter,beg);
			if(d1<d2)
				{vl=CC->nb;
				CC->list[vl]=nx;
				CC->nb=vl+1;
				exc[nx]=1;
				getf_find_rogc_todj(new_ter,&ter);
				}
			else if((d1>d2)&&(d2>eps_close))
				{vl=CC->nb;
				CC->list[vl]=nx;
				CC->nb=vl+1;
				exc[nx]=1;
				getf_find_rogc_todj(new_ter,&ter);
				}
			else
				break;
			}
		}
	}
return suc;
}



int qidk_conn_bowz(c_arc3D *ca,
int nb,hash_entry *CC,double eps_close)
{int i,j,*exc,n_cc,suc;
exc=(int *)malloc(nb*sizeof(int));
for(i=0;i<nb;i++)
	exc[i]=0;
n_cc=0;
for(j=0;j<2*nb;j++)
	{suc=mipt_curr_tucf(ca,nb,exc,&CC[n_cc],eps_close);
	n_cc++;
	if(suc==FAILURE)
		break;
	}
free(exc);
return n_cc;
}


int juwt_fork_hupt(c_arc3D *CA,int nb,int *excl,int *p_1,int *p_2)
{int i,j,p,q,suc;
double sml,dis,D;
point *A,*B;
A=(point *)malloc(2*sizeof(point));
B=(point *)malloc(2*sizeof(point));
suc=FAILURE;
sml=LARGE_NUMBER;
for(i=0;i<nb;i++)if(excl[i]==0)
for(j=0;j<i;j++)if(excl[j]==0)
	{getf_find_rogc_todj(CA[i].begn,&A[0]);
	getf_find_rogc_todj(CA[i].term,&A[1]);
	getf_find_rogc_todj(CA[j].begn,&B[0]);
	getf_find_rogc_todj(CA[j].term,&B[1]);
	D=LARGE_NUMBER;
	for(p=0;p<2;p++)
	for(q=0;q<2;q++)
		{dis=wodt_dist_gilq(A[p],B[q]);
		if(dis<D)
			D=dis;
		}
	if(D<sml)
		{*p_1=i;
		*p_2=j;
		sml=D;
		suc=SUCCESS;
		}
	}
free(B);
free(A);
return suc;
} 


int pejs_test_potn(c_arc3D *CA,int nb,
c_arc3D ca,int *id)
{int res,i,ts;
double eps_end=1.0e-4;
double eps_mid=1.0e-3;
double eps_rad=1.0e-3;
double eps_cent=1.0e-3;
*id=-1;
res=0;
for(i=0;i<nb;i++)
	{ts=fuqc_test_kuml(CA[i],ca,eps_end,eps_mid,eps_rad,eps_cent);
	if(ts==1)
		{*id=i;
		res=1;
		break;
		}
	}
return res;
}



int gikl_atom_vumh(double probe,c_arc3D *ca,
int n_ca,trmsrf *surf)
{int suc,k,i,*excl,p_1,p_2,nb,sk,ts_ex,id;
c_arc3D joint,C1,C2,*CA;
nb=n_ca;
excl=(int *)malloc(2*n_ca*sizeof(int));
CA=(c_arc3D *)malloc(2*n_ca*sizeof(c_arc3D));
for(i=0;i<n_ca;i++)
	poms_find_resk_lonb(ca[i],&CA[i]);
for(i=0;i<2*n_ca;i++)
	excl[i]=0;
k=0;
for(i=0;i<2*n_ca;i++)
	{sk=juwt_fork_hupt(CA,nb,excl,&p_1,&p_2);
	if(sk==FAILURE)
		break;
	if(sk==SUCCESS)
		{poms_find_resk_lonb(CA[p_1],&C1);
		poms_find_resk_lonb(CA[p_2],&C2);
		suc=bows_fork_wevc(probe,C1,C2,&surf[k],&joint);
		if(suc==SUCCESS)
			{excl[p_1]=1;
			excl[p_2]=1;
			ts_ex=pejs_test_potn(CA,nb,joint,&id);
			if(ts_ex==0)
				{poms_find_resk_lonb(joint,&CA[nb]);
				excl[nb]=0;
				nb++;
				}
			else
				excl[id]=1;
			k++;
			}
		}
	}
free(excl);
free(CA);
return k;
}



int lern_unat_moln(trmsrf *surf3,prat_main_blend G,
int p,c_arc3D *CA,int *idx)
{int att_alpha,att_beta,att_gamma;
int val,j,e,n1,n2,w,nb,sd1,sd2,sd;
val=G.N[p].val;
if(G.N[p].type_node==Q_TYPE)
	{fprintf(tmpout,"This applies only to T_TYPE node\n");
	exit(0);
	}
if(val==3)
	return 0;

att_alpha=0;
att_beta =0;
att_gamma=0;
for(j=0;j<val;j++)
	{e=G.N[p].inc_edge[j];
	
	n1=G.E[e].frvrt;
	sd1=G.E[e].side_1;
	
	n2=G.E[e].scvrt;
	sd2=G.E[e].side_2;
	
	if(n1==p)		sd=sd1;
	else if(n2==p)  sd=sd2;
	else
		{fprintf(tmpout,"Incident kt not connected to node\n");
		exit(0);
		}
	if(sd==ON_ALPHA)	att_alpha=1;
	if(sd==ON_BETA)		att_beta =1;
	if(sd==ON_GAMMA)	att_gamma=1;
	}

w=G.N[p].trim_idx;
nb=0;
if(att_alpha==0)
	{poms_find_resk_lonb(surf3[w].pc.alpha,&CA[nb]);
	idx[nb]=w;
	nb++;
	}
if(att_beta==0)
	{poms_find_resk_lonb(surf3[w].pc.beta,&CA[nb]);
	idx[nb]=w;
	nb++;
	}
if(att_gamma==0)
	{poms_find_resk_lonb(surf3[w].pc.gamma,&CA[nb]);
	idx[nb]=w;
	nb++;
	}
return nb;
}


