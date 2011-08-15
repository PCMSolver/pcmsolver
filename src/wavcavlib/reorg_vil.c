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
#include "pln_sph.h"
#include "coarsequad.h"



void huvc_find_defw_jutc(double *val,int N,int *pos)
{int i,j,k,tp,*ord;
double lrg;
ord=(int *)malloc(N*sizeof(int));
for(i=0;i<N;i++)
	ord[i]=i;
for(k=N;k>=1;k--)
	{lrg=-LARGE_NUMBER;
	for(i=0;i<k;i++)
		{if(val[ord[i]]>lrg)
			{j=i;
			lrg=val[ord[i]];
			}
		}
	tp=ord[k-1];   ord[k-1]=ord[j];   ord[j]=tp;
	}
for(i=0;i<N;i++)
	{tp=ord[i];
	pos[tp]=i;
	}
free(ord);
}



void tevh_dich_hegt(double *val,int N,int piv_id,double *val1,int *map1,
int *N1,double *val2,int *map2,int *N2)
{int i,k1,k2;
k1=0;	k2=0;
for(i=0;i<N;i++)
	{if(i!=piv_id)
		{if(val[i]<=val[piv_id])
			{val1[k1]=val[i];
			map1[k1]=i;
			k1++;
			}
		else
			{val2[k2]=val[i];
			map2[k2]=i;
			k2++;
			}
		}
	else
		{val1[k1]=val[i];
		map1[k1]=i;
		k1++;
		}
	}
*N1=k1;
*N2=k2;
}


int cojr_choo_suzq(double *val,int N,int *exc,int n)
{int res=-1,i,ts,dummy;
for(i=0;i<N;i++)
	{ts=gonl_arra_govj(exc,n,i,&dummy);
	if(ts==0)
		{res=i;
		break;
		}
	}
if(res==-1)
	res=0;
return res;
}


void vitp_dich_cefj(double *val,int N,double *val1,int *map1,
int *N1,double *val2,int *map2,int *N2)
{int n1,n2,piv_id,nb_trials=5,p,ind,piv_lrg=-1,*exc,n_exc;
double ratio,d1,d2,min_rat=0.6,lrg;
exc=(int *)malloc(nb_trials*sizeof(int));
ind=1;	lrg=0;	n_exc=0;
for(p=0;p<nb_trials;p++)
	{piv_id=cojr_choo_suzq(val,N,exc,n_exc);
	tevh_dich_hegt(val,N,piv_id,val1,map1,&n1,val2,map2,&n2);
	d1=(double)n1;
	d2=(double)n2;
	if(d1<d2)	ratio=d1/d2;
	else		ratio=d2/d1;
	if(ratio>min_rat)
		{*N1=n1;
		*N2=n2;
		ind=2;
		break;
		}
	if(ratio>lrg)
		{lrg=ratio;
		piv_lrg=piv_id;
		}
	exc[n_exc]=piv_id;
	n_exc++;
	}
free(exc);
if(ind==2)
	return;
if(ind==1)
	{if(piv_lrg==-1)
		piv_lrg=0;
	tevh_dich_hegt(val,N,piv_lrg,val1,map1,N1,val2,map2,N2);
	}
}



void levp_adap_nirz(double *val,int N,int len_short,int *pos)
{int N1,N2,*map1,*map2,*pos1,*pos2;
int i,s;
double *val1,*val2;
if(N<=len_short)
	huvc_find_defw_jutc(val,N,pos);
else
	{val1=(double *)malloc((N+1)*sizeof(double));
	val2 =(double *)malloc((N+1)*sizeof(double));
	map1 =(int *)malloc((N+1)*sizeof(int));
	map2 =(int *)malloc((N+1)*sizeof(int));
	pos1 =(int *)malloc((N+1)*sizeof(int));
	pos2 =(int *)malloc((N+1)*sizeof(int));
	vitp_dich_cefj(val,N,val1,map1,&N1,val2,map2,&N2);
	if((N1>=1)&&(N2>=1))
		{levp_adap_nirz(val1,N1,len_short,pos1);
		levp_adap_nirz(val2,N2,len_short,pos2);
		for(i=0;i<N1;i++)
			{s=pos1[i];
			pos[map1[i]]=s;
			}
		for(i=0;i<N2;i++)
			{s=pos2[i];
			pos[map2[i]]=s+N1;
			}
		}
	else
		huvc_find_defw_jutc(val,N,pos);
	free(val1);	free(val2);
	free(map1);	free(map2);
	free(pos1);	free(pos2);
	}
}









