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
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"



void vufn_find_cupl_jipw(c_curve *cc,int n,int N,mult_conn *mc)
{int i,j,k,l,m;
parm *p;
mc->nb_inner_polygons=n-1;
mc->nb_vr_outer=N*cc[0].N;
for(i=1;i<n;i++)
	mc->nb_vr_inner[i-1]=N*cc[i].N;

p=(parm *)malloc(N*sizeof(parm));
k=0;
for(i=0;i<n;i++)
	{m=cc[i].N;
	for(j=0;j<m;j++)
		{tegn_disc_likp(cc[i].nc[j],N,p);
		
		for(l=0;l<N;l++)
			{cunl_find_qedf_rewn(p[l],&mc->vertex[k]);
			k++;
			}
		}
	}
free(p);
mc->v_grs=k;
}


void welc_allo_dubg(int nin,int nvt,mult_conn *mc)
{mc->vertex=(parm *)malloc(nvt*sizeof(parm));
mc->zt   =(double *)malloc(nvt*sizeof(double));
mc->flag   =(int *)malloc(nvt*sizeof(int));
mc->mapglob=(int *)malloc(nvt*sizeof(int));
mc->nb_vr_inner=(int *)malloc(nin*sizeof(int));
}
 

void saqw_dest_kiqf(mult_conn *mc)
{free(mc->vertex);
free(mc->zt);
free(mc->flag);
free(mc->mapglob);
free(mc->nb_vr_inner);
}


void wazd_expo_qevt(char *filename,parm *V,int N)
{int i;
FILE *fp;
fp=fopen(filename,"w");
for(i=0;i<N;i++)
	fprintf(fp,"%f  %f\n",V[i].u,V[i].v);
fprintf(fp,"%f  %f\n",V[0].u,V[0].v);
fclose(fp);
}


void mors_expo_nugt(mult_conn mc)
{int N1,N2,i;
parm *temp;
N1=mc.nb_vr_outer;
temp=(parm *)malloc(N1*sizeof(parm));
for(i=0;i<N1;i++)
	cunl_find_qedf_rewn(mc.vertex[i],&temp[i]);
wazd_expo_qevt("MIOTY/out_vert.dat",temp,N1);
free(temp);

if(mc.nb_inner_polygons==1)
	{N2=mc.nb_vr_inner[0];
	temp=(parm *)malloc(N2*sizeof(parm));
	for(i=0;i<N2;i++)
		cunl_find_qedf_rewn(mc.vertex[N1+i],&temp[i]);
	wazd_expo_qevt("MIOTY/in_vert.dat",temp,N2);
	free(temp);
	}
}







