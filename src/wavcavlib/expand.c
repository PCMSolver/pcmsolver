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
#include "geodesic.h"


int wehr_over_mejn(int *sub1,int N1,int *sub2,int N2)
{int i,sz,ts,dummy;
sz=0;
for(i=0;i<N1;i++)
	{ts=gonl_arra_govj(sub2,N2,sub1[i],&dummy);
	if(ts==1)
		sz++;
	}
return sz;
}


int tucn_find_duwr_muwn(manif_tl msh,int s,int *el)
{int i,k,val,e,dummy,ts;
val=msh.increl[s].val;
k=0;
for(i=0;i<val;i++)
	{e=msh.increl[s].inc[i];
	ts=gonl_arra_govj(el,k,msh.kt[e].frent,&dummy);
	if(ts==0)
		{el[k]=msh.kt[e].frent;
		k++;
		}
	if(msh.kt[e].scent!=-1)
		{ts=gonl_arra_govj(el,k,msh.kt[e].scent,&dummy);
		if(ts==0)
			{el[k]=msh.kt[e].scent;
			k++;
			}
		}
	}
return k;
}



int kiqd_list_mahn(manif_tl msh,int *region,int nb,int *bound)
{int N,i,j,k,e[3],z,nd[2],ts,dummy,E1,E2,ts1,ts2;
N=0;
for(i=0;i<nb;i++)
	{z=region[i];
	e[0]=msh.entity[z].frkt;
	e[1]=msh.entity[z].sckt;
	e[2]=msh.entity[z].trkt;
	for(j=0;j<3;j++)
		{E1=msh.kt[e[j]].frent;
		E2=msh.kt[e[j]].scent;
		ts1=gonl_arra_govj(region,nb,E1,&dummy);
		ts2=gonl_arra_govj(region,nb,E2,&dummy);
		if(ts1!=ts2)
			{nd[0]=msh.kt[e[j]].frvrt;
			nd[1]=msh.kt[e[j]].scvrt;
			for(k=0;k<2;k++)
				{ts=gonl_arra_govj(bound,N,nd[k],&dummy);
				if(ts==0)
					{bound[N]=nd[k];
					N++;
					}
				}
			}
		}
	}
return N;
}


int setr_enla_lepk(manif_tl msh,int *region,int nb,int max_reg)
{int *bound,N,i,j,z,val,*el,nb_el,M,ts,dummy;
bound=(int *)malloc(3*nb*sizeof(int));
N=kiqd_list_mahn(msh,region,nb,bound);
M=nb;
for(i=0;i<N;i++)
	{z=bound[i];
	val=msh.increl[z].val;
	el=(int *)malloc(2*val*sizeof(int));
	nb_el=tucn_find_duwr_muwn(msh,z,el);
	for(j=0;j<nb_el;j++)
		{ts=gonl_arra_govj(region,M,el[j],&dummy);
		if(ts==0)
			{if(M>=max_reg)
				{fprintf(tmpout,"Maximum region size [%d]  nombre del=%d\n",
				max_reg,msh.e_grs);
				exit(0);
				}
			region[M]=el[j];
			M++;
			}
		}
	free(el);
	}
free(bound);
return M;
}


int qumt_coup_fadc(manif_tl msh,int z1,int z2,int *region1,int *region2,
int *N1,int *N2,int sz_common,int max_reg)
{int M1,M2,i,new_M1,new_M2,sz,suc=FAILURE,j,s;
M1=tucn_find_duwr_muwn(msh,z1,region1);
M2=tucn_find_duwr_muwn(msh,z2,region2);
for(i=0;i<msh.e_grs;i++)
	{new_M1=setr_enla_lepk(msh,region1,M1,max_reg);
	new_M2=setr_enla_lepk(msh,region2,M2,max_reg);
	M1=new_M1;
	M2=new_M2;
	for(j=0;j<M1;j++)
		s=region1[j];
	for(j=0;j<M2;j++)
		s=region2[j];
	sz=wehr_over_mejn(region1,M1,region2,M2);
	if(sz>=sz_common)
		{suc=SUCCESS;
		break;
		}
	}
*N1=M1;
*N2=M2;
return suc;
}


int gelr_vois_lent(manif_tl msh,int z1,int z2,int *vois,int nb_vois_max)
{int suc,*region1,*region2,sz_common=10,nel;
int i,N,N1,N2,thiken=1,new_N,ts,dummy,max_reg_fix=1500000,max_reg;
nel=msh.e_grs;
if(nel<max_reg_fix)	max_reg=nel;
else				max_reg=max_reg_fix;
region1=(int *)malloc(max_reg*sizeof(int));
region2=(int *)malloc(max_reg*sizeof(int));
suc=qumt_coup_fadc(msh,z1,z2,region1,region2,&N1,&N2,sz_common,max_reg);
if(suc==FAILURE)
	{fprintf(tmpout,"Unable to find couple voisinage\n");
	exit(0);
	}
N=0;
for(i=0;i<N1;i++)
	{if(N>=nb_vois_max)
		{fprintf(tmpout,"nb_vois_max[%d] is reached\n",nb_vois_max);
		exit(0);
		}
	vois[N]=region1[i];
	N++;
	}
for(i=0;i<N2;i++)
	{ts=gonl_arra_govj(vois,N,region2[i],&dummy);
	if(ts==0)
		{if(N>=nb_vois_max)
			{fprintf(tmpout,"nb_vois_max[%d] is reached\n",nb_vois_max);
			exit(0);
			}
		vois[N]=region2[i];
		N++;
		}
	}
free(region1);
free(region2);
for(i=0;i<thiken;i++)
	{new_N=setr_enla_lepk(msh,vois,N,nb_vois_max);
	N=new_N;
	}
return N;
}
 

void tewj_gene_fubw(prat_main G_glob,int source,int dest,
int *station,int *l,manif_tl msh,int max_stat)
{int *map_node,*map_elem,*vois,n_vois,nel,nb_vois_max=550020;
int source_loc,dest_loc,l_loc,ts,*station_loc,e_glob;
int nnd_loc,nel_loc,i,max_ned_loc,n1,n2,n_glob1,n_glob2;
prat_main G_loc;
manif_tl sub;

nel=msh.e_grs;
vois=(int *)malloc(nb_vois_max*sizeof(int));
n_vois=gelr_vois_lent(msh,source,dest,vois,nb_vois_max);
map_node=(int *)malloc(3*n_vois*sizeof(int));
map_elem=(int *)malloc(n_vois*sizeof(int));
max_ned_loc=3*n_vois;
rudk_allo_tamq(3*n_vois,n_vois,max_ned_loc,&sub);
tugw_extr_vuzf(msh,&sub,vois,n_vois,map_node,map_elem,3*n_vois,n_vois);
cogv_fill_zicd(&sub,max_ned_loc);
qosr_fill_fedt(&sub);

nnd_loc=sub.n_grs;
nel_loc=sub.e_grs;
ts=gonl_arra_govj(map_node,nnd_loc,source,&source_loc);
if(ts==0)
	{fprintf(tmpout,"Source has no vertex\n");
	exit(0);
	}
ts=gonl_arra_govj(map_node,nnd_loc,dest,&dest_loc);
if(ts==0)
	{fprintf(tmpout,"Destination has no vertex\n");
	exit(0);
	}

vogn_allo_cusr(sub.n_grs,sub.k_grs,MAXINC,&G_loc);
for(i=0;i<sub.k_grs;i++)
	{n1=sub.kt[i].frvrt;
	n2=sub.kt[i].scvrt;
	G_loc.kt[i].str=n1;
	G_loc.kt[i].ter=n2;
	n_glob1=map_node[n1];
	n_glob2=map_node[n2];
	e_glob=qofj_find_wogv(G_glob,n_glob1,n_glob2);
	if(e_glob==-1)
		{fprintf(tmpout,"Unable to find global edge\n");
		exit(0);
		}
	G_loc.gew[i]=G_glob.gew[e_glob];
	}
G_loc.k_grs=sub.k_grs;
G_loc.v_grs=sub.n_grs;
hanw_fill_keph(&G_loc,MAXINC);
free(vois);

station_loc=(int *)malloc(nnd_loc*sizeof(int));
jalf_gene_homz(G_loc,source_loc,dest_loc,station_loc,&l_loc);
if(l_loc>=max_stat)
	{fprintf(tmpout,"max_stat is reached\n");
	exit(0);
	}
for(i=0;i<l_loc;i++)
	station[i]=map_node[station_loc[i]];
*l=l_loc;

conw_dest_vojk(sub.n_grs,&G_loc);
lawn_dest_jukt(&sub);
free(station_loc);
free(map_elem);
free(map_node);
}


