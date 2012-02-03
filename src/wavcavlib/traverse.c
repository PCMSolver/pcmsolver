
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
#include "splinemol.h"
#include "pln_sph.h"
#include "sas.h"
#include "geodesic.h"



void fuch_fuse_nuqt(manif_tl mshin,manif_tl *mshout,int *map_node)
{int nnd,nel,nb_pcw,*ls_pcw,i,j;
int s,ind,*map,n1,n2,n3,*MU,ts;
double eps=1.0e-6;

nnd=mshin.n_grs;
nel=mshin.e_grs;
map=(int *)malloc(nnd*sizeof(int));
for(i=0;i<nnd;i++)
	map[i]=i;

ls_pcw=(int *)malloc(nnd*sizeof(int));
nb_pcw=1;
ls_pcw[0]=0;
map[0]=0;
map_node[0]=0;
for(i=1;i<nnd;i++)
	{ind=1;
	for(j=0;j<nb_pcw;j++)
		{s=ls_pcw[j];
		ts=gect_tole_husn(mshin.knot[s],mshin.knot[i],eps);
		if(ts==1)
			{ind=2;
			map[i]=s;
			map_node[i]=s;
			break;
			}
		}
	if(ind==1)
		{s=nb_pcw;
		ls_pcw[s]=i;
		map_node[i]=s;
		map[i]=i;
		nb_pcw++;
		}
	}



MU=(int *)malloc(nnd*sizeof(int));
mshout->n_grs=nb_pcw;
for(i=0;i<nb_pcw;i++)
	{s=ls_pcw[i];
	MU[s]=i;
	mshout->knot[i].absi=mshin.knot[s].absi;
	mshout->knot[i].ordo=mshin.knot[s].ordo;
	mshout->knot[i].cote=mshin.knot[s].cote;
	}
free(ls_pcw);

mshout->e_grs=nel;
for(i=0;i<nel;i++)
	{n1=mshin.entity[i].frvrt;	s=map[n1];
	mshout->entity[i].frvrt=MU[s];
	n2=mshin.entity[i].scvrt;	s=map[n2];
	mshout->entity[i].scvrt=MU[s];
	n3=mshin.entity[i].thvrt;	s=map[n3];
	mshout->entity[i].thvrt=MU[s];
	}
free(MU);
free(map);
}



void nowj_fuse_cogs(manif_tl *msh,int *map_node)
{int nnd,nel,i,n1,n2,n3;
manif_tl mshout;

nnd=msh->n_grs;
nel=msh->e_grs;
mshout.knot   =(point *)malloc(nnd*sizeof(point));
mshout.entity=(telolf *)malloc(nel*sizeof(telolf));
fuch_fuse_nuqt(*msh,&mshout,map_node);

nnd=mshout.n_grs;
for(i=0;i<nnd;i++)
	{msh->knot[i].absi=mshout.knot[i].absi;
	msh->knot[i].ordo=mshout.knot[i].ordo;
	msh->knot[i].cote=mshout.knot[i].cote;
	}
nel=mshout.e_grs;
for(i=0;i<nel;i++)
	{n1=mshout.entity[i].frvrt;
	n2=mshout.entity[i].scvrt;
	n3=mshout.entity[i].thvrt;
	msh->entity[i].frvrt=n1;
	msh->entity[i].scvrt=n2;
	msh->entity[i].thvrt=n3;
	}
msh->n_grs=nnd;
msh->e_grs=nel;
free(mshout.knot);
free(mshout.entity);
}



int rumq_same_defg(int n1,int n2,int n3,int m1,int m2,int m3)
{int nd[3],md[3],idx,i,j,res;
nd[0]=n1;	nd[1]=n2;	nd[2]=n3;
md[0]=m1;	md[1]=m2;	md[2]=m3;
res=1;  idx=0;
for(i=0;i<3;i++)
	{idx=0;
	for(j=0;j<3;j++)
		{if(nd[i]==md[j])
			{idx=1;
			break;
			}
		}
	if(idx==0)
		{res=0;
		break;
		}
	}
return res;
}


int pazj_test_zucn(telolf t1,telolf t2)
{int n1,n2,n3,res,m1,m2,m3;
n1=t1.frvrt;
n2=t1.scvrt;
n3=t1.thvrt;
m1=t2.frvrt;
m2=t2.scvrt;
m3=t2.thvrt;
res=rumq_same_defg(n1,n2,n3,m1,m2,m3);
return res;
}


void tevm_fuse_nocr(manif_tl *msh,rgb_lk *col,sphere *S)
{int nel,i,j,k,ts,ind;
telolf *temp;
rgb_lk *col_temp;
sphere *S_temp;
nel=msh->e_grs;
temp=(telolf *)malloc(nel*sizeof(telolf));
col_temp=(rgb_lk *)malloc(nel*sizeof(rgb_lk));
S_temp=(sphere *)malloc(nel*sizeof(sphere));
k=0;
for(i=0;i<nel;i++)
	{ind=1;
	for(j=0;j<k;j++)
		{ts=pazj_test_zucn(msh->entity[i],temp[j]);
		if(ts==1)
			{ind=2;
			break;
			}
		}
	if(ind==1)
		{temp[k].frvrt=msh->entity[i].frvrt;
		temp[k].scvrt =msh->entity[i].scvrt;
		temp[k].thvrt =msh->entity[i].thvrt;
		col_temp[k].red  =col[i].red;
		col_temp[k].green=col[i].green;
		col_temp[k].blue =col[i].blue;
		neqg_find_lodr_bogm(S[i],&S_temp[k]);
		k++;
		}
	}
for(i=0;i<k;i++)
	{msh->entity[i].frvrt=temp[i].frvrt;
	msh->entity[i].scvrt=temp[i].scvrt;
	msh->entity[i].thvrt=temp[i].thvrt;
	col[i].red  =col_temp[i].red;
	col[i].green=col_temp[i].green;
	col[i].blue =col_temp[i].blue;
	neqg_find_lodr_bogm(S_temp[i],&S[i]);
	}
msh->e_grs=k;


free(temp);
free(col_temp);
free(S_temp);
}


void leqd_inte_wags(manif_tl *msh,rgb_lk *col,sphere *S,int t,
manif_tl loc,int max_nnd,int max_nel)
{int nnd,nel,NND,NEL,i,k,*map,n1,n2,n3;
nnd=loc.n_grs;
nel=loc.e_grs;
NND=msh->n_grs;
NEL=msh->e_grs;

map=(int *)malloc(nnd*sizeof(int));
k=NND;
for(i=0;i<nnd;i++)
	{if(k>=max_nnd)
		{fprintf(tmpout,"Allocated memory [%d] for nodes is not enough\n",max_nnd);
		exit(0);
		}
	getf_find_rogc_todj(loc.knot[i],&msh->knot[k]);
	map[i]=k;
	k++;
	}
msh->n_grs=k;


if(t>=max_nel)
	{fprintf(tmpout,"Allocated memory for elements is not enough\n");
	exit(0);
	}
n1=loc.entity[0].frvrt;	msh->entity[t].frvrt=map[n1];
n2=loc.entity[0].scvrt;	msh->entity[t].scvrt=map[n2];
n3=loc.entity[0].thvrt;	msh->entity[t].thvrt=map[n3];

k=NEL;
for(i=1;i<nel;i++)
	{if(k>=max_nel)
		{fprintf(tmpout,"Allocated memory for elements is not enough\n");
		exit(0);
		}
	n1=loc.entity[i].frvrt;	msh->entity[k].frvrt=map[n1];
	n2=loc.entity[i].scvrt;	msh->entity[k].scvrt=map[n2];
	n3=loc.entity[i].thvrt;	msh->entity[k].thvrt=map[n3];
	col[k].red  =col[t].red;
	col[k].green=col[t].green;
	col[k].blue =col[t].blue;
	neqg_find_lodr_bogm(S[t],&S[k]);
	k++;
	}
msh->e_grs=k;
free(map);
}


int lekf_refi_jech(manif_tl *msh,rgb_lk *col,sphere *S,
trav_tri *T,int nb,int max_nnd,int max_nel)
{int i,j,N,r,n[4],*type_s,*type_t,suc,SUC;
point *A,*s,*t;
manif_tl loc;
A=(point *)malloc(3*sizeof(point));
for(i=0;i<nb;i++)
	{r=T[i].el_idx;
	n[1]=msh->entity[r].frvrt;
	n[2]=msh->entity[r].scvrt;
	n[3]=msh->entity[r].thvrt;
	for(j=1;j<=3;j++)
		getf_find_rogc_todj(msh->knot[n[j]],&A[j-1]);
	N=T[i].nb_trav;
	s=(point *)malloc(N*sizeof(point));
	t=(point *)malloc(N*sizeof(point));
	type_s=(int *)malloc(N*sizeof(int));
	type_t=(int *)malloc(N*sizeof(int));
	for(j=0;j<N;j++)
		{getf_find_rogc_todj(T[i].v_str[j],&s[j]);
		getf_find_rogc_todj(T[i].v_ter[j],&t[j]);
		type_s[j]=T[i].type_str[j];
		type_t[j]=T[i].type_trm[j];
		}
	
	loc.knot=(point *)malloc(10*(N+1)*sizeof(point));
	loc.entity=(telolf *)malloc(10*(N+1)*sizeof(telolf));
	suc=fils_loca_poth(A,s,t,type_s,type_t,N,&loc);
	if(suc==FAILURE)
		SUC=FAILURE;
	else
	  {leqd_inte_wags(msh,col,S,r,loc,max_nnd,max_nel);
	    SUC=SUCCESS;
	  }
	free(loc.knot);
	free(loc.entity);
	free(s);		free(t);
	free(type_s);	free(type_t);
	if(SUC==FAILURE)
		break;
	}
free(A);
return SUC;
}



int pojw_list_rasq(manif_tl msh,float_curve *FC,
int N,trav_tri *T)
{int nb,i,j,M,cas,ts,t,idx,*tau,nel,s;
nel=msh.e_grs;
tau=(int *)malloc(nel*sizeof(int));
nb=0;
for(i=0;i<N;i++)
	{M=FC[i].st_grs;
	for(j=0;j<M-1;j++)
		{cas=FC[i].cs[j];
		if((2<=cas)&&(cas<=4))
			{t=FC[i].trv[j];
			ts=gonl_arra_govj(tau,nb,t,&idx);
			if(ts==0)
				{tau[nb]=t;
				T[nb].el_idx=t;
				T[nb].nb_trav=1;
				getf_find_rogc_todj(FC[i].stn[j],&T[nb].v_str[0]);
				getf_find_rogc_todj(FC[i].stn[j+1],&T[nb].v_ter[0]);
				if(cas==2)
					{T[nb].type_str[0]=APEX_POINT;
					T[nb].type_trm[0]=EDGE_POINT;
					}
				if(cas==3)
					{T[nb].type_str[0]=EDGE_POINT;
					T[nb].type_trm[0]=EDGE_POINT;
					}
				if(cas==4)
					{T[nb].type_str[0]=EDGE_POINT;
					T[nb].type_trm[0]=APEX_POINT;
					}
				nb++;
				}
			else
				{s=T[idx].nb_trav;
				getf_find_rogc_todj(FC[i].stn[j],&T[idx].v_str[s]);
				getf_find_rogc_todj(FC[i].stn[j+1],&T[idx].v_ter[s]);
					if(cas==2)
					{T[idx].type_str[s]=APEX_POINT;
					T[idx].type_trm[s]=EDGE_POINT;
					}
				if(cas==3)
					{T[idx].type_str[s]=EDGE_POINT;
					T[idx].type_trm[s]=EDGE_POINT;
					}
				if(cas==4)
					{T[idx].type_str[s]=EDGE_POINT;
					T[idx].type_trm[s]=APEX_POINT;
					}
				T[idx].nb_trav=s+1;
				}
			}
		}
	}
free(tau);
return nb;
}


void fukl_disp_talf(float_curve fc)
{int i;
fprintf(tmpout,"Number of stations=%d\n",fc.st_grs);
for(i=0;i<fc.st_grs;i++)
	fprintf(tmpout,"station[%d]=[%f,%f,%f]\n",i,fc.stn[i].absi,
	fc.stn[i].ordo,fc.stn[i].cote);
for(i=0;i<fc.st_grs;i++)
	fprintf(tmpout,"i=%d  nd_idx=%d  ed_idx=%d  [%f,%f,%f]\n",
	i,fc.nd_idx[i],fc.kt_idx[i],fc.stn[i].absi,
	fc.stn[i].ordo,fc.stn[i].cote);
fprintf(tmpout,"\n");
for(i=0;i<fc.st_grs-1;i++)
	fprintf(tmpout,"i=%d  case=%d  trav=%d\n",i,fc.cs[i],fc.trv[i]);
fprintf(tmpout,"----------------\n");
}


void rekf_find_kujn_rukg(float_curve fc_in,float_curve *fc_out)
{int N,i;
N=fc_in.st_grs;
for(i=0;i<N;i++)
	{getf_find_rogc_todj(fc_in.stn[i],&fc_out->stn[i]);
	fc_out->nd_idx[i] =fc_in.nd_idx[i];
	fc_out->kt_idx[i] =fc_in.kt_idx[i];
	fc_out->nd_stat[i]=fc_in.nd_stat[i];
	}
for(i=0;i<N-1;i++)
	{fc_out->cs[i]      =fc_in.cs[i];
	fc_out->trv[i]=fc_in.trv[i];
	fc_out->kt_stat[i]  =fc_in.kt_stat[i];
	}
fc_out->st_grs=N;
}



int danp_comp_laqz(manif_tl msh,float_curve *fc,int k)
{int e,n1,n2,ind,q,i,suc=FAILURE;
double eps=1.0e-6,dis;
point *A;
float_curve temp;
A=(point *)malloc(2*sizeof(point));
if(fc->trv[k]==fc->trv[k+1])
	{if((fc->cs[k]==2)&&(fc->cs[k+1]==4))
		{e=fc->kt_idx[k+1];
		n1=msh.kt[e].frvrt;
		n2=msh.kt[e].scvrt;
		getf_find_rogc_todj(msh.knot[n1],&A[0]);
		getf_find_rogc_todj(msh.knot[n2],&A[1]);
		ind=1;
		for(i=0;i<2;i++)
			{dis=wodt_dist_gilq(A[i],fc->stn[k+2]);
			if(dis<eps)
				{ind=2;
				break;
				}
			}
		if(ind==2)
			{vewk_allo_jovk(fc->st_grs-1,&temp);
			q=0;
			for(i=0;i<fc->st_grs;i++)if(i!=k+1)
				{getf_find_rogc_todj(fc->stn[i],&temp.stn[q]);
				temp.nd_idx[q]=fc->nd_idx[i];
				temp.kt_idx[q]=fc->kt_idx[i];
				q++;
				}
			for(i=0;i<k;i++)
				{temp.cs[i]=fc->cs[i];
				temp.trv[i]=fc->trv[i];
				}
			temp.cs[k]=1;
			temp.trv[k]=fc->trv[k+1];
			for(i=k+2;i<fc->st_grs-1;i++)
				{temp.cs[i-1]=fc->cs[i];
				temp.trv[i-1]=fc->trv[i];
				}
			temp.st_grs=fc->st_grs-1;
			rekf_find_kujn_rukg(temp,fc);
			lohm_dest_nosr(&temp);
			suc=SUCCESS;
			}
		}
	}
free(A);
return suc;
}


int lajt_comp_dipr(manif_tl msh,float_curve *fc)
{int k,suc,SUC=FAILURE;
for(k=0;k<fc->st_grs-2;k++)
	{suc=danp_comp_laqz(msh,fc,k);
	if(suc==SUCCESS)
		{SUC=SUCCESS;
		break;
		}
	}
return SUC;
}


void zekl_redu_govr(manif_tl msh,float_curve *fc)
{int suc,i,N;
N=fc->st_grs;
for(i=0;i<N;i++)
	{suc=lajt_comp_dipr(msh,fc);
	if(suc==FAILURE)
		break;
	}
}



 
