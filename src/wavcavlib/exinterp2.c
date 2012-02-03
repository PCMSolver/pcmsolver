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



int nawd_unat_pifr(included_sides *inc,trmsrf *surf2,
prat_main_blend G,int p,c_arc3D *CA,int *idx)
{int att_alpha,att_beta,att_gamma,att_delta,nb_inc;
int val,j,e,n1,n2,w,nb,sd1,sd2,sd,on_sd[3];
int not_inc[2],k,ts,dummy;
int qw[4]={ON_ALPHA,ON_BETA,ON_GAMMA,ON_DELTA};
val=G.N[p].val;
if(G.N[p].type_node==T_TYPE)
	{fprintf(tmpout,"This applies only to Q_TYPE node\n");
	exit(0);
	}
if(val==2)
	return 0;
not_inc[0]=inc[p].side_1;
not_inc[1]=inc[p].side_2;
k=1;
for(j=0;j<4;j++)
	{ts=gonl_arra_govj(not_inc,2,qw[j],&dummy);
	if(ts==0)
		{on_sd[k]=qw[j];
		k++;
		}
	}

att_alpha=0;	att_beta =0;
att_gamma=0;	att_delta=0;
nb_inc=0;
if((on_sd[1]==ON_ALPHA)||(on_sd[2]==ON_ALPHA))
	{att_alpha=1;
	nb_inc++;
	}
if((on_sd[1]==ON_BETA)||(on_sd[2]==ON_BETA))
	{att_beta=1;
	nb_inc++;
	}
if((on_sd[1]==ON_GAMMA)||(on_sd[2]==ON_GAMMA))
	{att_gamma=1;
	nb_inc++;
	}
if((on_sd[1]==ON_DELTA)||(on_sd[2]==ON_DELTA))
	{att_delta=1;
	nb_inc++;
	}
if(nb_inc!=2)
	fprintf(tmpout,"WARNING: not 2 sides are incident upon atoms\n");
for(j=0;j<val;j++)
	{e=G.N[p].inc_edge[j];
	
	n1 =G.E[e].frvrt;
	sd1=G.E[e].side_1;
	
	n2 =G.E[e].scvrt;
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
	if(sd==ON_DELTA)	att_delta=1;
	}

w=G.N[p].trim_idx;
nb=0;
if(att_alpha==0)
	{poms_find_resk_lonb(surf2[w].pt.alpha,&CA[nb]);
	idx[nb]=w;
	nb++;
	}
if(att_beta==0)
	{poms_find_resk_lonb(surf2[w].pt.beta,&CA[nb]);
	idx[nb]=w;
	nb++;
	}
if(att_gamma==0)
	{poms_find_resk_lonb(surf2[w].pt.gamma,&CA[nb]);
	idx[nb]=w;
	nb++;
	}
if(att_delta==0)
	{poms_find_resk_lonb(surf2[w].pt.delta,&CA[nb]);
	idx[nb]=w;
	nb++;
	}
return nb;
}



int celf_unat_nozp(included_sides *inc,trmsrf *surf2,
trmsrf *surf3,prat_main_blend G,int p,
c_arc3D *CA,int *idx)
{int nb;
if(G.N[p].type_node==T_TYPE)
	nb=lern_unat_moln(surf3,G,p,CA,idx);
else if(G.N[p].type_node==Q_TYPE)
	nb=nawd_unat_pifr(inc,surf2,G,p,CA,idx);
else
	{fprintf(tmpout,"p=%d    type=%d\n",p,G.N[p].type_node);
	fprintf(tmpout,"p-th G-node is neither T-type nor Q-type\n");
	exit(0);
	}
return nb;
}


int dejc_unat_hibn(included_sides *inc,
trmsrf *surf2,trmsrf *surf3,
prat_main_blend G,c_arc3D *CA,int *idx)
{int nnd,i,j,n_loc,nb,*idx_loc;
c_arc3D *ca_loc;
nnd=G.n_grs;
ca_loc=(c_arc3D *)malloc(4*sizeof(c_arc3D));
idx_loc=(int *)malloc(4*sizeof(int));
nb=0;
for(i=0;i<nnd;i++)
	{n_loc=celf_unat_nozp(inc,surf2,surf3,G,i,ca_loc,idx_loc);
	for(j=0;j<n_loc;j++)
		{poms_find_resk_lonb(ca_loc[j],&CA[nb]);
		CA[nb].c_cir=0;
		idx[nb]=idx_loc[j];
		nb++;
		}
	}
free(idx_loc);
free(ca_loc);
return nb;
}


void caqw_disp_pegr(prat_main_blend G)
{int nnd,ned,i,n1,n2;
nnd=G.n_grs;
fprintf(tmpout,"Nb incidence nodes=%d\n",nnd);
for(i=0;i<nnd;i++)
	{if(G.N[i].type_node==Q_TYPE)
		fprintf(tmpout,"n[%d]  Q trim_idx=%d   valence=%d\n",i,G.N[i].trim_idx,G.N[i].val);
	if(G.N[i].type_node==T_TYPE)
		fprintf(tmpout,"n[%d]  T trim_idx=%d   valence=%d\n",i,G.N[i].trim_idx,G.N[i].val);
	}
ned=G.k_grs;
fprintf(tmpout,"Nb incidence edges=%d\n",ned);
for(i=0;i<ned;i++)
	{n1=G.E[i].frvrt;
	n2=G.E[i].scvrt;
	fprintf(tmpout,"e[%d]   [n1,n2]=[%d,%d]  [side_1,side_2]=",i,n1,n2);
	if(G.E[i].side_1==ON_ALPHA)
		fprintf(tmpout,"[ON_ALPHA,");
	if(G.E[i].side_1==ON_BETA)
		fprintf(tmpout,"[ON_BETA,");
	if(G.E[i].side_1==ON_GAMMA)
		fprintf(tmpout,"[ON_GAMMA,");
	if(G.E[i].side_1==ON_DELTA)
		fprintf(tmpout,"[ON_DELTA,");
	
	if(G.E[i].side_2==ON_ALPHA)
		fprintf(tmpout,"ON_ALPHA]\n");
	if(G.E[i].side_2==ON_BETA)
		fprintf(tmpout,"ON_BETA]\n");
	if(G.E[i].side_2==ON_GAMMA)
		fprintf(tmpout,"ON_GAMMA]\n");
	if(G.E[i].side_2==ON_DELTA)
		fprintf(tmpout,"ON_DELTA]\n");
	}
}



int qumw_atom_picm(atom *S,double probe,trmsrf *surf1,
int nb_surf1,trmsrf *surf2,int nb_surf2,trmsrf *surf3,
int nb_surf3,set_arcs *SA,blend_cpx BC,double eps_close,
int max_surf3,hash_entry *gap,int *n_gap,int max_gap,int max_gap_val)
{int i,j,z,nb,nb_cc,M,p,n_loc,nnd,n_gp;
int N,*idx,*idx_temp,n_worst,n_new;
prop_ccurve pcc;
c_arc3D *ca,*C,*CA;
trmsrf *loc;
included_sides *inc;
hash_entry *CC;
prat_main_blend G;

nnd=0;
for(i=0;i<nb_surf2;i++)
	{if(surf2[i].type==4)
		nnd++;
	if(surf2[i].type==6)
		{fprintf(tmpout,"We dont treat type six here\n");
		exit(0);
		}
	}
nnd=nnd+nb_surf3;
N=BC.bt_grs;
inc=(included_sides *)malloc(N*sizeof(included_sides));
G.N=(node_blend *)malloc(nnd*sizeof(node_blend));
G.E=(kt_blend *)malloc(3*nnd*sizeof(kt_blend));
det_graph_bl_simple(S,surf2,nb_surf2,surf3,nb_surf3,BC,&G,inc);


n_worst=4*nb_surf2+3*nb_surf3;
CA=(c_arc3D *)malloc(n_worst*sizeof(c_arc3D));
idx_temp=(int *)malloc(n_worst*sizeof(int));
nb=dejc_unat_hibn(inc,surf2,surf3,G,CA,idx_temp);
ca=(c_arc3D *)malloc(nb*sizeof(c_arc3D));
idx=(int *)malloc(nb*sizeof(int));
for(i=0;i<nb;i++)
	{poms_find_resk_lonb(CA[i],&ca[i]);
	idx[i]=idx_temp[i];
	}
free(idx_temp);	free(CA);
free(G.N);		free(G.E);
free(inc);		free(idx);


CC=(hash_entry *)malloc(nb*sizeof(hash_entry));
for(i=0;i<nb;i++)
	CC[i].list=(int *)malloc(nb*sizeof(int));
nb_cc=qidk_conn_bowz(ca,nb,CC,eps_close);
n_new=nb_surf3; 
n_gp=0;
for(p=0;p<nb_cc;p++)
	{M=CC[p].nb;
	C=(c_arc3D *)malloc(M*sizeof(c_arc3D));
	for(j=0;j<M;j++)
		{z=CC[p].list[j];
		poms_find_resk_lonb(ca[z],&C[j]);
		}
	qezj_find_tukd_wesg(&pcc);
	loc=(trmsrf *)malloc(M*sizeof(trmsrf));
	for(i=0;i<M;i++)
		{homd_allo_tevf(pcc,&loc[i].cc);
		surf3[i].pc.sitn=0;
		homd_allo_tevf(pcc,&loc[i].pc.cc);
		}
	
	n_loc=gikl_atom_vumh(probe,C,M,loc);
	free(C);
	for(i=0;i<n_loc;i++)
		{if(n_new>=max_surf3)
			{fprintf(tmpout,"max_surf3 is reached\n");
			exit(0);
			}
		jofd_find_mikn_gehj(loc[i].pc,&surf3[n_new].pc);		
		surf3[n_new].nb_inner=0;
		kotg_find_wuhk_kemt(loc[i].cc,&surf3[n_new].cc);
		surf3[n_new].type    =5;
		surf3[n_new].boundary=1;
		surf3[n_new].nb_inner=0;
		surf3[n_new].pc.sitn =0;
		if(i>=max_gap_val)
			{fprintf(tmpout,"max_gap_val=%d is reached\n",max_gap_val);
			exit(0);
			}
		if(n_gp>=max_gap)
			{fprintf(tmpout,"max_gap=%d is reached\n",max_gap);
			exit(0);
			}
		gap[n_gp].list[i]=n_new;
		n_new++;
		}
	gap[n_gp].nb=n_loc;
	n_gp++;
	for(i=0;i<M;i++)
		{wosn_dest_jomw(pcc,&loc[i].cc);
		wosn_dest_jomw(pcc,&loc[i].pc.cc);
		}
	free(loc);
	}
*n_gap=n_gp;

for(i=0;i<nb;i++)
	free(CC[i].list);
free(CC);
free(ca);
return n_new;
}

 

