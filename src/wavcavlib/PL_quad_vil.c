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
#include "geodesic.h"
#include "smooth.h"
#include "pln_sph.h"
#include "sas.h"
#include "coarsequad.h"


int bi_valence_obsol(int ob,fajor_sion3D quad)
{int i,nel,n1,n2,n3,n4,val,res=1;
nel=quad.e_grs;
val=0;
for(i=0;i<nel;i++)
	{n1=quad.elem[i].frvrt;
	n2=quad.elem[i].scvrt;
	n3=quad.elem[i].thvrt;
	n4=quad.elem[i].ftvrt;
	if((n1==ob)||(n2==ob)||(n3==ob)||(n4==ob))
		{val++;
		if(val>=3)
			{res=0;
			break;
			}
		}
	}
if(val!=2)
	res=0;
return res;
}


int verif_quad_closed(fajor_sion3D quad)
{int ned,i,e1,e2,ind,res=1;
ned=quad.k_grs;
ind=1;
for(i=0;i<ned;i++)
	{e1=quad.kt[i].frent;
	e2=quad.kt[i].scent;
	if(e2==-1)
		{fprintf(tmpout,"kt=%d  [e1,e2]=[%d,%d]\n",i,e1,e2);
		res=0;
		ind=2;
		break;
		}
	}
return res;
}



int nusg_merg_ribq(fajor_sion3D quad,int *p,int *q,int *obs)
{int ned,i,E,F,ts,ob,nb,ts_bi;
int tr1,tr2,tr3,tr4,dummy;
ned=quad.k_grs;
nb=0;
for(i=0;i<ned;i++)
	{E=quad.kt[i].frent;
	F=quad.kt[i].scent;
	if(F==-1)
		fprintf(tmpout,"WARNING: quadrangulation is not closed\n");
	if(F!=-1)
		{ts=cojs_shar_sejm(quad,E,F,&ob);
		if(ts==1)
			{tr1=gonl_arra_govj(p,nb,E,&dummy);
			tr2=gonl_arra_govj(p,nb,F,&dummy);
			tr3=gonl_arra_govj(q,nb,E,&dummy);
			tr4=gonl_arra_govj(q,nb,F,&dummy);
			if((tr1==0)&&(tr2==0)&&(tr3==0)&&(tr4==0))
				{p[nb]=E;
				q[nb]=F;
				ts_bi=bi_valence_obsol(ob,quad);
				if(ts_bi==1)
					{fprintf(tmpout,"obolete node=%d   couple quads=[%d,%d]\n",ob,E,F);
					obs[nb]=ob;
					nb++;
					}
				}
			}
		}
	}
return nb;
}



int rech_astr_qonk(fajor_sion3D quad,int p,int q,int obs)
{int res=-1,nd[4],md[4],i,ts,dummy;
nd[0]=quad.elem[p].frvrt;
nd[1]=quad.elem[p].scvrt;
nd[2]=quad.elem[p].thvrt;
nd[3]=quad.elem[p].ftvrt;
md[0]=quad.elem[q].frvrt;
md[1]=quad.elem[q].scvrt;
md[2]=quad.elem[q].thvrt;
md[3]=quad.elem[q].ftvrt;
for(i=0;i<4;i++)
	{ts=gonl_arra_govj(nd,4,md[i],&dummy);
	if((ts==0)&&(md[i]!=obs))
		{res=md[i];
		break;
		}
	}
if(res==-1)
	{jans_disp_nudj(quad,p);
	jans_disp_nudj(quad,q);
	fprintf(tmpout,"obs=%d\n",obs);
	fprintf(tmpout,"Unable to find astride node\n");
	}
return res;
}


void vunc_merg_levm(fajor_sion3D *quad,int p,int q,
int ob,int *obs_node,int *n_nd,int *obs_elem,int *n_el)
{int as_nd,k;
as_nd=rech_astr_qonk(*quad,p,q,ob);
if(quad->elem[p].frvrt==ob)	quad->elem[p].frvrt=as_nd;
if(quad->elem[p].scvrt==ob)	quad->elem[p].scvrt=as_nd;
if(quad->elem[p].thvrt==ob)	quad->elem[p].thvrt=as_nd;
if(quad->elem[p].ftvrt==ob)	quad->elem[p].ftvrt=as_nd;
k=*n_nd;
obs_node[k]=ob;
*n_nd=k+1;
k=*n_el;
obs_elem[k]=q;
*n_el=k+1;
}


void durs_disc_vohq(fajor_sion3D *quad,int *obs,int M)
{int nnd,i,k,ts,dummy,*map,n1,n2,n3,n4;
point *temp;
nnd=quad->n_grs;
temp=(point *)malloc(nnd*sizeof(point));
map=(int *)malloc(nnd*sizeof(int));
k=0;
for(i=0;i<nnd;i++)
	{ts=gonl_arra_govj(obs,M,i,&dummy);
	if(ts==0)
		{getf_find_rogc_todj(quad->knot[i],&temp[k]);
		map[i]=k;
		k++;
		}
	}
for(i=0;i<quad->e_grs;i++)
	{n1=quad->elem[i].frvrt;	quad->elem[i].frvrt=map[n1];
	n2=quad->elem[i].scvrt;	quad->elem[i].scvrt=map[n2];
	n3=quad->elem[i].thvrt;	quad->elem[i].thvrt=map[n3];
	n4=quad->elem[i].ftvrt;	quad->elem[i].ftvrt=map[n4];
	}
for(i=0;i<k;i++)
	getf_find_rogc_todj(temp[i],&quad->knot[i]);
quad->n_grs=k;
free(temp);
free(map);
}


void rosv_disc_vinl(fajor_sion3D *quad,int *el,int M)
{int i,k,nel,ts,dummy;
efajor *temp;
nel=quad->e_grs;
temp=(efajor *)malloc(nel*sizeof(efajor));
k=0;
for(i=0;i<nel;i++)
	{ts=gonl_arra_govj(el,M,i,&dummy);
	if(ts==0)
		{temp[k].frvrt=quad->elem[i].frvrt;
		temp[k].scvrt=quad->elem[i].scvrt;
		temp[k].thvrt=quad->elem[i].thvrt;
		temp[k].ftvrt=quad->elem[i].ftvrt;
		k++;
		}
	}

quad->e_grs=nel-M;
for(i=0;i<nel-M;i++)
	{quad->elem[i].frvrt=temp[i].frvrt;
	quad->elem[i].scvrt=temp[i].scvrt;
	quad->elem[i].thvrt=temp[i].thvrt;
	quad->elem[i].ftvrt=temp[i].ftvrt;
	}
free(temp);
}


int seqg_find_gihc_secb(fajor_sion3D *quad,int max_ned,
int *force_term)
{int *p,*q,*obs,ned,nb,i,suc=FAILURE;
int *obs_node,*obs_elem,n_nd,n_el,ts;
*force_term=0;
ts=verif_quad_closed(*quad);
if(ts==0)
	{fprintf(tmpout,"quadrangulation is not closed\n");
	*force_term=1;
	return FAILURE;
	}
ned=quad->k_grs;
p  =(int *)malloc(ned*sizeof(int));
q  =(int *)malloc(ned*sizeof(int));
obs=(int *)malloc(ned*sizeof(int));
nb =nusg_merg_ribq(*quad,p,q,obs);
fprintf(tmpout,"Number of mergeables=%d\n",nb);

if(nb>=1)
	{suc=SUCCESS;
	obs_node=(int *)malloc(nb*sizeof(int));
	obs_elem=(int *)malloc(nb*sizeof(int));
	n_nd=0;
	n_el=0;
	for(i=0;i<nb;i++)
		vunc_merg_levm(quad,p[i],q[i],obs[i],obs_node,&n_nd,obs_elem,&n_el);
	durs_disc_vohq(quad,obs_node,n_nd);
	rosv_disc_vinl(quad,obs_elem,n_el);
	free(obs_node);
	free(obs_elem);
	teqr_fill_regm(quad,max_ned);
	}
free(q);
free(obs);
return suc;
}


void detl_find_faqs_cakq(fajor_sion3D *quad,
int max_ned,int *forc_term)
{int suc,i,nnd,f_trm;
*forc_term=0;
nnd=quad->n_grs;
for(i=0;i<nnd;i++)
	{suc=seqg_find_gihc_secb(quad,max_ned,&f_trm);
	if(f_trm==1)
		{*forc_term=1;
		return;
		}
	if(suc==FAILURE)
		break;
	verif_quad_closed(*quad);
	}
}



int jold_dete_nojq(fajor_sion3D quad,int z,double scale,int *nd_a,int *nd_b)
{int n[4],res;
double D1,D2;
n[0]=quad.elem[z].frvrt;
n[1]=quad.elem[z].scvrt;
n[2]=quad.elem[z].thvrt;
n[3]=quad.elem[z].ftvrt;
D1=wodt_dist_gilq(quad.knot[n[0]],quad.knot[n[2]]);
D2=wodt_dist_gilq(quad.knot[n[1]],quad.knot[n[3]]);
res=0;
if(D2>scale*D1)
	{res=1;
	*nd_a=n[0];
	*nd_b=n[2];
	}
else if(D1>scale*D2)
	{res=1;
	*nd_a=n[1];
	*nd_b=n[3];
	}
return res;
}


void tefc_disc_jazq(fajor_sion3D *quad,
int z,int nd_a,int nd_b)
{int nel,i,n1,n2,n3,n4,*el,*obs;
point X;
X.absi=0.5*(quad->knot[nd_a].absi+quad->knot[nd_b].absi);
X.ordo=0.5*(quad->knot[nd_a].ordo+quad->knot[nd_b].ordo);
X.cote=0.5*(quad->knot[nd_a].cote+quad->knot[nd_b].cote);
getf_find_rogc_todj(X,&quad->knot[nd_a]);
nel=quad->e_grs;
for(i=0;i<nel;i++)
	{n1=quad->elem[i].frvrt;	if(n1==nd_b)  quad->elem[i].frvrt=nd_a;
	n2=quad->elem[i].scvrt;	if(n2==nd_b)  quad->elem[i].scvrt=nd_a;
	n3=quad->elem[i].thvrt;	if(n3==nd_b)  quad->elem[i].thvrt=nd_a;
	n4=quad->elem[i].ftvrt;	if(n4==nd_b)  quad->elem[i].ftvrt=nd_a;
	}
el=(int *)malloc(sizeof(int));
el[0]=z;
rosv_disc_vinl(quad,el,1);
free(el);
obs=(int *)malloc(sizeof(int));
obs[0]=nd_b;
durs_disc_vohq(quad,obs,1);
free(obs);
}


void ciql_remo_cekq(fajor_sion3D *quad,
double scale,int max_ned)
{int z,nd_a,nd_b,ts,ind,nb_trials=3,i,change=0;

for(i=0;i<nb_trials;i++)
	{ind=1;
	for(z=0;z<quad->e_grs;z++)
		{ts=jold_dete_nojq(*quad,z,scale,&nd_a,&nd_b);
		if(ts==1)
			{tefc_disc_jazq(quad,z,nd_a,nd_b);
			ind=2;
			change=1;
			}
		}
	if(ind==1)
		break;
	}
if(change==1)
	teqr_fill_regm(quad,max_ned);

}

