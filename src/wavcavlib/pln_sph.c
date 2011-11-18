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
#include <math.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"



int bosr_plan_revg(point G,vect3D nrm,sphere S,
point *omega,double *rad)
{double dis,sp,R;
vect3D temp;
bofp_form_nukv(S.zent,G,&temp);
sp=rocv_scal_toqc(temp,nrm);
dis=fabs(sp);
R=S.rad;
if(dis>R)
	return 0;

omega->absi=S.zent.absi+sp*nrm.absi;
omega->ordo=S.zent.ordo+sp*nrm.ordo;
omega->cote=S.zent.cote+sp*nrm.cote;
*rad=sqrt(R*R-dis*dis);
return 1;
} 


int peqj_test_vesf(circle3D c,sphere S)
{int res,Z;
double rad,D;
point omega;
res=0;
Z=bosr_plan_revg(c.zent,c.nrml,S,&omega,&rad);
if(Z==1)
	{D=wodt_dist_gilq(omega,c.zent);
	if(D+c.rad<rad)
		res=1;
	}
return res;
}


void berd_disp_derh(sphere S)
{fprintf(tmpout,"center=[%f,%f,%f]\n",S.zent.absi,S.zent.ordo,S.zent.cote);
fprintf(tmpout,"radius=%f\n",S.rad);
}



int pavz_circ_kuts(circle3D c,int N,sphere *S,int nb_S)
{int q,i,j,n_loc,ind,res,ts,*Z;
double step,t,phi,theta,*rad,d1,d2,d3;
c_arc3D *c_loc;
point *omega,mid,*P;

ind=1;
for(i=0;i<nb_S;i++)
	{ts=peqj_test_vesf(c,S[i]);
	if(ts==1)
		{ind=2;
		break;
		}
	}
if(ind==2)
	return 1;

vewr_sphe_ruhd(c.nrml.absi,c.nrml.ordo,
c.nrml.cote,&phi,&theta);
step=1.0/(double)N;
P=(point *)malloc(N*sizeof(point));
for(q=0;q<N;q++)
	{t=step*(double)q;
	mesl_punc_guvf(c,t,&P[q]);
	}
c_loc=(c_arc3D *)malloc(N*sizeof(c_arc3D));
dept_spli_piwg(c,phi,theta,P,N,c_loc);
n_loc=N;
free(P);


omega=(point *)malloc(nb_S*sizeof(point));
rad=(double *)malloc(nb_S*sizeof(double));
Z=(int *)malloc(nb_S*sizeof(int));
for(i=0;i<nb_S;i++)
	Z[i]=bosr_plan_revg(c.zent,c.nrml,S[i],&omega[i],&rad[i]);

res=1;
for(j=0;j<n_loc;j++)
	{puwj_midp_curq(c_loc[j],&mid);	
	ind=1;
	for(i=0;i<nb_S;i++)if(Z[i]==1)
		{d1=wodt_dist_gilq(omega[i],mid);
		if(d1<rad[i])
			{d2=wodt_dist_gilq(omega[i],c_loc[j].begn);
			if(d2<rad[i])
				{d3=wodt_dist_gilq(omega[i],c_loc[j].term);
				if(d3<rad[i])
					{ind=2;
					break;
					}
				}
			}
		}
	if(ind==1)
		{res=0;
		break;
		}
	}


free(rad);
free(omega);
free(c_loc);
free(Z);
return res;
}



