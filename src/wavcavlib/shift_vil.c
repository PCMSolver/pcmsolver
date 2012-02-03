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
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "smooth.h"


int dikr_shif_qigf(prat_main GR,manif_tl msh,float_curve fc,
int N1,int N2,int *list,int max_list)
{int nb,i,j,ts,dummy,*nd_arr,N,t,nd[3],n_stat;
n_stat=fc.st_grs;
nd_arr=(int *)malloc(3*n_stat*sizeof(int));
N=0;
for(i=0;i<n_stat-1;i++)
	{t=fc.trv[i];
	nd[0]=msh.entity[t].frvrt;
	nd[1]=msh.entity[t].scvrt;
	nd[2]=msh.entity[t].thvrt;
	for(j=0;j<3;j++)
		{ts=gonl_arra_govj(nd_arr,N,nd[j],&dummy);
		if(ts==0)
			{nd_arr[N]=nd[j];
			N++;
			}
		}
	}
nd[0]=N1;
nd[1]=N2;
for(j=0;j<2;j++)
	{ts=gonl_arra_govj(nd_arr,N,nd[j],&dummy);
	if(ts==0)
		{nd_arr[N]=nd[j];
		N++;
		}
	}

nb=porw_vici_nejt(GR,msh,nd_arr,N,list,max_list);
free(nd_arr);
return nb;
}



int qocj_shif_vukc(prat_main GR,manif_tl msh,int nb_inter,float_curve fc,
int N1,int N2,point *omega,int *type,int *underl,int max_stat_smooth)
{int nb,*list,R,nel;
int max_list_aux=300,max_list;

nel=msh.e_grs;
if(nel<max_list_aux)  max_list=nel;
else				  max_list=max_list_aux;
list=(int *)malloc(max_list*sizeof(int));
R=dikr_shif_qigf(GR,msh,fc,N1,N2,list,max_list);

nb=hubs_smoo_citg(GR,msh,N1,N2,list,R,nb_inter,omega,
type,underl,max_stat_smooth);
free(list);
return nb;
}


void zumc_shif_rumc(prat_main GR,manif_tl msh,int nb_inter,
float_curve fc,int N1,int N2,float_curve *fc_new,int nb_new)
{int nb_stat=1,*type,*underl,max_stat=20000;
point *omega;
omega=(point *)malloc(max_stat*sizeof(point));
type=(int *)malloc(max_stat*sizeof(int));
underl=(int *)malloc(max_stat*sizeof(int));
nb_stat=qocj_shif_vukc(GR,msh,nb_inter,fc,
N1,N2,omega,type,underl,max_stat);

if(nb_stat>=nb_new)
	{fprintf(tmpout,"nb_new is reached\n");
	exit(0);
	}
with_conv_davf(msh,GR,omega,type,underl,nb_stat,fc_new);
free(omega);
free(type);
free(underl);
}


int mesc_cond_wonp(fajor_sion3D *QUAD,int *corresp,hash_entry *H,
manif_tl msh,vect3D *nrm,prat_main GR,int nb_inter,int z,int w,float_curve *F,int *new_w)
{int i,*ngb,val,n1,n2,suc=FAILURE,nb_new=100;
int e,w_new,q,N1,N2,str,trm,nb_st,E;
point X;
double cur,sft,lrg;
float_curve fc_new;

val=GR.dgr[w];
ngb=(int *)malloc(val*sizeof(int));
for(i=0;i<val;i++)
	{e=GR.incd[w][i];
	n1=GR.kt[e].ter;
	n2=GR.kt[e].str;
	if(n1==w)	ngb[i]=n2;
	if(n2==w)	ngb[i]=n1;
	}

cur=lord_qual_weqf(z,*QUAD,corresp,nrm,F,H);
lrg=-LARGE_NUMBER;
for(i=0;i<val;i++)
	{w_new=ngb[i];
	getf_find_rogc_todj(msh.knot[w_new],&X);
	sft=rohl_qual_zifh(z,X,nrm[w_new],*QUAD,corresp,msh,nrm,F,H);
	if(sft>lrg)
		{lrg=sft;
		q=i;
		}
	}
w_new=ngb[q];
getf_find_rogc_todj(msh.knot[w_new],&X);
free(ngb);

if((lrg>cur)&&(lrg>0.0))
	{getf_find_rogc_todj(X,&QUAD->knot[z]);
	for(i=0;i<H[z].nb;i++)
		{E=H[z].list[i];
		nb_st=F[E].st_grs;
		str=F[E].nd_idx[0];
		trm=F[E].nd_idx[nb_st-1];
		if(str==w)
			{N1=w_new;
			N2=trm;
			}
		else if(trm==w)
			{N1=w_new;
			N2=str;
			}
		else
			{fprintf(tmpout,"Float curve must have w\n");
			exit(0);
			}
		vewk_allo_jovk(nb_new,&fc_new);
		zumc_shif_rumc(GR,msh,nb_inter,F[E],N1,N2,&fc_new,nb_new);
		rekf_find_kujn_rukg(fc_new,&F[E]);
		lohm_dest_nosr(&fc_new);
		}
	*new_w=w_new;
	suc=SUCCESS;
	}
return suc;
}



void hopf_impr_zugw(fajor_sion3D *QUAD,manif_tl msh,vect3D *nrm,
int *corresp,prat_main GR,float_curve *F,int nb_times)
{int max_loc_inc=20,i,nnd_q,z,w,new_w,suc;
int nb_inter=5,p;
hash_entry *H;
nnd_q=QUAD->n_grs;
H=(hash_entry *)malloc(nnd_q*sizeof(hash_entry));
for(i=0;i<nnd_q;i++)
	H[i].list=(int *)malloc(max_loc_inc*sizeof(int));
piln_inci_lung(*QUAD,H,max_loc_inc);

for(p=0;p<nb_times;p++)
for(z=0;z<nnd_q;z++)
	{fprintf(tmpout,"SHIFT QUAD VERTEX[%d/%d]  traversal=%d\n",z,nnd_q-1,p);
	
	w=corresp[z];
	suc=mesc_cond_wonp(QUAD,corresp,H,msh,nrm,GR,nb_inter,z,w,F,&new_w);
	if(suc==SUCCESS)
		corresp[z]=new_w;
	
	}
for(i=0;i<nnd_q;i++)
	free(H[i].list);
free(H);
}


void jans_disp_nudj(fajor_sion3D quad,int z)
{int e[4],i;
fprintf(tmpout,"quadrilateral[%d]\n",z);
fprintf(tmpout,"nodes=[%d,%d,%d,%d]\n",quad.elem[z].frvrt,quad.elem[z].scvrt,
quad.elem[z].thvrt,quad.elem[z].ftvrt);
e[0]=quad.elem[z].frkt;
e[1]=quad.elem[z].sckt;
e[2]=quad.elem[z].trkt;
e[3]=quad.elem[z].ftkt;
fprintf(tmpout,"edges=[%d,%d,%d,%d]\n",e[0],e[1],e[2],e[3]);
for(i=0;i<4;i++)
	fprintf(tmpout,"loc edge[%d]=%d  nodes=[%d,%d]  elems=[%d,%d]\n",i,e[i],
	quad.kt[e[i]].frvrt,quad.kt[e[i]].scvrt,quad.kt[e[i]].frent,quad.kt[e[i]].scent);
fprintf(tmpout,"-------------\n");
}
 


