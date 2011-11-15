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


double merc_leng_jugb(point *omega,int N)
{int i;
double L;
L=0.0;
for(i=0;i<N-1;i++)
	L=L+wodt_dist_gilq(omega[i],omega[i+1]);
return L;
}


double fujt_leng_nebd(manif_tl msh,int *station,int N)
{int i,s;
double res;
point *omega;
omega=(point *)malloc(N*sizeof(point));
for(i=0;i<N;i++)
	{s=station[i];
	omega[i].absi=msh.knot[s].absi;
	omega[i].ordo=msh.knot[s].ordo;
	omega[i].cote=msh.knot[s].cote;
	}
res=merc_leng_jugb(omega,N);
free(omega);
return res;
}



int hubs_smoo_citg(prat_main GR,manif_tl msh,int N1,int N2,int *list,
int R,int N,point *omega,int *type,int *underl,int max_stat_smooth)
{int *under,n1,n2,i,nb_mu,*mu,L,*stat,s,im,n_vert,val;
double len;
point *coord;
prat_main g_loc;
coord=(point *)malloc(3*N*R*sizeof(point));
under=(int *)malloc(3*N*R*sizeof(int));
mu=(int *)malloc(3*R*sizeof(int));

n_vert=(3+3*N)*R;
vogn_allo_cusr(n_vert,n_vert*MAXINC,MAXINC,&g_loc);
nb_mu=pahb_vici_sinj(N,GR,list,R,msh,&g_loc,coord,under,mu);
hanw_fill_keph(&g_loc,MAXINC);

val=0;
for(i=0;i<g_loc.v_grs;i++)
if(g_loc.dgr[i]>val)
	val=g_loc.dgr[i];
if(val>=MAXINC)
	{fprintf(tmpout,"MAXINC is reached\n");
	exit(0);
	}

n1=-1;		n2=-1;
for(i=0;i<nb_mu;i++)
	{if(mu[i]==N1)	n1=i;
	if(mu[i]==N2)	n2=i;
	if((n1!=-1)&&(n2!=-1))
		break;
	}


stat=(int *)malloc(g_loc.v_grs*sizeof(int));
jalf_gene_homz(g_loc,n1,n2,stat,&L);
for(i=0;i<L;i++)
	{s=stat[i];
	if(s<nb_mu)
		{if(i>=max_stat_smooth)
			{fprintf(tmpout,"max_stat_smooth is had\n");
			exit(0);
			}
		type[i]=2;
		im=mu[s];
		omega[i].absi=msh.knot[im].absi;
		omega[i].ordo=msh.knot[im].ordo;
		omega[i].cote=msh.knot[im].cote;
		underl[i]=im;
		}
	else
		{if(i>=max_stat_smooth)
			{fprintf(tmpout,"max_stat_smooth is had\n");
			exit(0);
			}
		type[i]=1;
		omega[i].absi=coord[s-nb_mu].absi;
		omega[i].ordo=coord[s-nb_mu].ordo;
		omega[i].cote=coord[s-nb_mu].cote;
		underl[i]=under[s-nb_mu];
		}
	}
len=merc_leng_jugb(omega,L);
free(stat);		free(mu);
free(coord);	free(under);
conw_dest_vojk(n_vert,&g_loc);
return L;
}



int kusv_star_wuqg(int nd,prat_main GR,manif_tl msh,int *ngb)
{int nb,val,i,j,e,el[2];
nb=0;
val=GR.dgr[nd];
for(i=0;i<val;i++)
	{e=GR.incd[nd][i];
	el[0]=msh.kt[e].frent;
	el[1]=msh.kt[e].scent;
	for(j=0;j<2;j++)
	if(el[j]!=-1)
		{ngb[nb]=el[j];
		nb++;
		}
	}
return nb;
}



int porw_vici_nejt(prat_main GR,manif_tl msh,
int *station,int L,int *list,int max_list)
{int R,i,j,k,nb,*ngb,nd,ts,dummy;
ngb=(int *)malloc(msh.n_grs*sizeof(int));
k=0;
for(i=0;i<L;i++)
	{nd=station[i];
	nb=kusv_star_wuqg(nd,GR,msh,ngb);
	for(j=0;j<nb;j++)
		{ts=gonl_arra_govj(list,k,ngb[j],&dummy);
		if(ts==0)
			{if(k>=max_list)
				{fprintf(tmpout,"max_list=%d is reached\n",max_list);
				exit(0);
				}
			list[k]=ngb[j];
			k++;
			}
		}
	}
R=k;
free(ngb);
return R;
}



int sowt_expa_dacw(prat_main GR,manif_tl msh,int *list,
int L,int max_list)
{int N,i,j,l,s,n[3],*ngb,ts,dummy,nb;
N=L;
ngb=(int *)malloc(msh.n_grs*sizeof(int));
for(i=0;i<L;i++)
	{s=list[i];
	n[0]=msh.entity[s].frvrt;
	n[1]=msh.entity[s].scvrt;
	n[2]=msh.entity[s].thvrt;
	for(j=0;j<3;j++)
		{nb=kusv_star_wuqg(n[j],GR,msh,ngb);
		for(l=0;l<nb;l++)
			{ts=gonl_arra_govj(list,N,ngb[l],&dummy);
			if(ts==0)
				{if(N>=max_list)
					{fprintf(tmpout,"max_list=%d is had\n",max_list);
					exit(0);
					}
				list[N]=ngb[l];
				N++;
				}
			}
		}
	}
free(ngb);
return N;
}



int leqg_path_joqw(prat_main GR,manif_tl msh,int *station,
int L,int level,int *list,int max_list)
{int nb,m,i,new_m;
m=porw_vici_nejt(GR,msh,station,L,list,max_list);
for(i=0;i<=level;i++)
	{new_m=sowt_expa_dacw(GR,msh,list,m,max_list);
	if(new_m==m)		
		break;
	m=new_m;
	}
nb=m;
return nb; 
}



int folc_gene_kost(int N1,int N2,int level,int N,prat_main GR,
manif_tl msh,point *omega,int *type,int *underl,int max_stat_smooth)
{int *list,nnd,nel,R,nb,nb_stat,*station,i,s;
int max_list_aux=100000,max_list,max_stat,max_stat_aux=20000;
double len;
nnd=msh.n_grs;
nel=msh.e_grs;
if(nnd<max_stat_aux)  max_stat=nnd;
else				  max_stat=max_list_aux;
station=(int *)malloc(max_stat*sizeof(int));
tewj_gene_fubw(GR,N1,N2,station,&nb_stat,msh,max_stat);
for(i=0;i<nb_stat;i++)
	s=station[i];
len=fujt_leng_nebd(msh,station,nb_stat);
if(nel<max_list_aux)  max_list=nel;
else				  max_list=max_list_aux;
list=(int *)malloc(max_list*sizeof(int));
R=leqg_path_joqw(GR,msh,station,nb_stat,level,list,max_list);
free(station);
nb=hubs_smoo_citg(GR,msh,N1,N2,list,R,N,omega,type,
underl,max_stat_smooth);
free(list);
return nb;
}
 


