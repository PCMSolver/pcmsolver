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
#include "sas.h"
#include "smooth.h"


void piln_inci_lung(fajor_sion3D QUAD,
hash_entry *H,int max_loc_inc)
{int i,j,nnd,ned,nd[2],val;
nnd=QUAD.n_grs;
ned=QUAD.k_grs;
for(i=0;i<nnd;i++)
	H[i].nb=0;
for(j=0;j<ned;j++)
	{nd[0]=QUAD.kt[j].frvrt;
	nd[1]=QUAD.kt[j].scvrt;
	for(i=0;i<2;i++)
		{val=H[nd[i]].nb;
		if(val>=max_loc_inc)
			{fprintf(tmpout,"max_loc_inc is had\n");
			exit(0);
			}
		H[nd[i]].list[val]=j;
		H[nd[i]].nb=val+1;
		}
	}
}


void kugr_best_mogh(point *sam,int n_sam,point *X)
{int i;
double *len,t=0.5;
if(n_sam==1)
	getf_find_rogc_todj(sam[0],X);
else
	{len=(double *)malloc(n_sam*sizeof(double));
	len[0]=0.0;
	for(i=1;i<n_sam;i++)
		len[i]=len[i-1]+wodt_dist_gilq(sam[i-1],sam[i]);
	vefm_eval_bilt(t,sam,n_sam,len,X);
	free(len);
	}
}


void vuwr_best_lerq(point A,point B,manif_tl msh,
int *cand,int n,point *X)
{int i,z,n_sam=3;
double sml,d1,d2,d;
point *sam,temp;
sml=LARGE_NUMBER;
for(i=0;i<n;i++)
	{z=cand[i];
	d1=wodt_dist_gilq(A,msh.knot[z]);
	d2=wodt_dist_gilq(B,msh.knot[z]);
	d=d1+d2;
	if(d<sml)
		{sml=d;
		getf_find_rogc_todj(msh.knot[z],&temp);
		}
	}
sam=(point *)malloc(n_sam*sizeof(point));
getf_find_rogc_todj(A,&sam[0]);
getf_find_rogc_todj(temp,&sam[1]);
getf_find_rogc_todj(B,&sam[2]);
kugr_best_mogh(sam,n_sam,X);
free(sam);
}



void gupf_best_tawd(point A,point B,manif_tl msh,
float_curve F,point *X)
{int *cand,n_trav,i,j;
int n,t,nd[3],ts,dummy;
n_trav=F.st_grs-1;
cand=(int *)malloc(3*n_trav*sizeof(int));
n=0;
for(i=0;i<n_trav;i++)
	{t=F.trv[i];
	nd[0]=msh.entity[t].frvrt;
	nd[1]=msh.entity[t].scvrt;
	nd[2]=msh.entity[t].thvrt;
	for(j=0;j<3;j++)
		{ts=gonl_arra_govj(cand,n,nd[j],&dummy);
		if(ts==0)
			{cand[n]=nd[j];
			n++;
			}
		}
	}
vuwr_best_lerq(A,B,msh,cand,n,X);
free(cand);

}


void nopr_best_jiwt(float_curve F,point *X)
{int n_sam;
n_sam=F.st_grs;
kugr_best_mogh(F.stn,n_sam,X);
}


int weqf_edge_rikn(fajor_sion3D QUAD,int w,int n1,int n2)
{int *e,res=-1,i,nd_a,nd_b;
e=(int *)malloc(4*sizeof(int));
e[0]=QUAD.elem[w].frkt;
e[1]=QUAD.elem[w].sckt;
e[2]=QUAD.elem[w].trkt;
e[3]=QUAD.elem[w].ftkt;
for(i=0;i<4;i++)
	{nd_a=QUAD.kt[e[i]].frvrt;
	nd_b=QUAD.kt[e[i]].scvrt;
	if(((n1==nd_a)&&(n2==nd_b))||((n2==nd_a)&&(n1==nd_b)))
		{res=e[i];
		break;
		}
	}
free(e);
return res;
}


int lqs=0;



double zejc_loca_form(manif_tl msh,vect3D *nrm,float_curve *FC,
fajor_sion3D QUAD,int *corresp,int z,point X,vect3D nrm_X,int w)
{int nd[4],i,cr,nx,e,q;
double Q;
point *cor,*mid,A,B;
vect3D *N;

cor=(point *)malloc(4*sizeof(point));
N=(vect3D *)malloc(4*sizeof(vect3D));
nd[0]=QUAD.elem[w].frvrt;
nd[1]=QUAD.elem[w].scvrt;
nd[2]=QUAD.elem[w].thvrt;
nd[3]=QUAD.elem[w].ftvrt;
for(i=0;i<4;i++)
	{if(nd[i]==z)
		{getf_find_rogc_todj(X,&cor[i]);
		getf_find_rogc_todj(nrm_X,&N[i]);
		}
	else
		{getf_find_rogc_todj(QUAD.knot[nd[i]],&cor[i]);
		q=corresp[nd[i]];
		getf_find_rogc_todj(nrm[q],&N[i]);
		}
	}

mid=(point *)malloc(4*sizeof(point));
for(i=0;i<4;i++)
	{cr=nd[i];
	if(i==3)
		nx=nd[0];
	else
		nx=nd[i+1];
	e=weqf_edge_rikn(QUAD,w,cr,nx);
	if(e==-1)
		{fprintf(tmpout,"2-Unable to find local edge\n");
		exit(0);
		}
	if((cr!=z)&&(nx!=z))
		nopr_best_jiwt(FC[e],&mid[i]);
	else
		{getf_find_rogc_todj(X,&A);
		if(z==cr)
			getf_find_rogc_todj(QUAD.knot[nx],&B);
		else if(z==nx)
			getf_find_rogc_todj(QUAD.knot[cr],&B);
		gupf_best_tawd(A,B,msh,FC[e],&mid[i]);
		}
	}

Q=lurq_qual_wotl(cor,mid,N,0);
free(cor);
free(mid);
free(N);
return Q;
}



double lijf_loca_necf(float_curve *FC,fajor_sion3D QUAD,
int *corresp,vect3D *nrm,int w)
{int nd[4],i,cr,nx,e,z;
double Q;
point *cor,*mid;
vect3D *N;

cor=(point *)malloc(4*sizeof(point));
N=(vect3D *)malloc(4*sizeof(vect3D));
nd[0]=QUAD.elem[w].frvrt;
nd[1]=QUAD.elem[w].scvrt;
nd[2]=QUAD.elem[w].thvrt;
nd[3]=QUAD.elem[w].ftvrt;
for(i=0;i<4;i++)
	{getf_find_rogc_todj(QUAD.knot[nd[i]],&cor[i]);
	z=corresp[nd[i]];
	getf_find_rogc_todj(nrm[z],&N[i]);
	}

mid=(point *)malloc(4*sizeof(point));
for(i=0;i<4;i++)
	{cr=nd[i];
	if(i==3)
		nx=nd[0];
	else
		nx=nd[i+1];
	e=weqf_edge_rikn(QUAD,w,cr,nx);
	if(e==-1)
		{fprintf(tmpout,"[cr,nx]=[%d,%d]\n",cr,nx);
		fprintf(tmpout,"3-Unable to find local edge\n");
		exit(0);
		}
	nopr_best_jiwt(FC[e],&mid[i]);
	}

Q=lurq_qual_wotl(cor,mid,N,0);
free(N);
free(cor);
free(mid);
return Q;
}



double rohl_qual_zifh(int z,point X,vect3D nrm_X,fajor_sion3D QUAD,
int *corresp,manif_tl msh,vect3D *nrm,float_curve *FC,hash_entry *H)
{int *el,val,i,j,k,E,e[2];
int ts,dummy,n_inc,w;
double qual,Q;

val=H[z].nb;
el=(int *)malloc(2*val*sizeof(int));
k=0;
for(i=0;i<val;i++)
	{E=H[z].list[i];
	e[0]=QUAD.kt[E].frent;
	e[1]=QUAD.kt[E].scent;
	for(j=0;j<2;j++)if(e[j]!=-1)
		{ts=gonl_arra_govj(el,k,e[j],&dummy);
		if(ts==0)
			{el[k]=e[j];
			k++;
			}
		}
	}
n_inc=k;

qual=LARGE_NUMBER;
for(i=0;i<n_inc;i++)
	{w=el[i];
	Q=zejc_loca_form(msh,nrm,FC,QUAD,corresp,z,X,nrm_X,w);
	if(Q<=qual)
		qual=Q;
	}
free(el);
return qual;
}



double lord_qual_weqf(int z,fajor_sion3D QUAD,int *corresp,
vect3D *nrm,float_curve *FC,hash_entry *H)
{int *el,val,i,j,k,E,e[2];
int ts,dummy,n_inc,w;
double qual,Q;

val=H[z].nb;
el=(int *)malloc(2*val*sizeof(int));
k=0;
for(i=0;i<val;i++)
	{E=H[z].list[i];
	e[0]=QUAD.kt[E].frent;
	e[1]=QUAD.kt[E].scent;
	for(j=0;j<2;j++)
		{ts=gonl_arra_govj(el,k,e[j],&dummy);
		if(ts==0)
			{el[k]=e[j];
			k++;
			}
		}
	}
n_inc=k;

qual=LARGE_NUMBER;
for(i=0;i<n_inc;i++)
	{w=el[i];
	Q=lijf_loca_necf(FC,QUAD,corresp,nrm,w);
	if(Q<=qual)
		qual=Q;
	}
free(el);
return qual;
}



