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
#include <math.h>
#include <stdlib.h>
#include <stdlib.h>
#include "cavity.h"
#include "splinemol.h"
#include "pln_sph.h"
#include "coarsequad.h"


void jasg_rear_cupg(int ned,double *err,int *list)
{int i,len_short=20,*pos,s;
pos=(int *)malloc(ned*sizeof(int));
levp_adap_nirz(err,ned,len_short,pos);
for(i=0;i<ned;i++)
	{s=pos[i];
	list[s]=i;
	}
free(pos);
}



void folp_fuse_ruqk(manif_tl *sub,int m1,int m2,int e,point P)
{int nnd,nel,i,k,el1,el2,n1,n2,n3,*map;
double srf;
telolf *temp_el;
point *temp_nd;

nnd=sub->n_grs;
temp_nd=(point *)malloc(nnd*sizeof(point));
map=(int *)malloc(nnd*sizeof(int));
k=0;
for(i=0;i<nnd;i++)if((i!=m1)&&(i!=m2))
	{temp_nd[k].absi=sub->knot[i].absi;
	temp_nd[k].ordo=sub->knot[i].ordo;
	temp_nd[k].cote=sub->knot[i].cote;
	map[i]=k;
	k++;
	}
temp_nd[k].absi=P.absi;
temp_nd[k].ordo=P.ordo;
temp_nd[k].cote=P.cote;
map[m1]=k;
map[m2]=k;
k++;
sub->n_grs=k;

el1=sub->kt[e].frent;
el2=sub->kt[e].scent;
nel=sub->e_grs;
temp_el=(telolf *)malloc(nel*sizeof(telolf));
k=0;
for(i=0;i<nel;i++)if((i!=el1)&&(i!=el2))
	{n1=sub->entity[i].frvrt;	 temp_el[k].frvrt=map[n1];
	n2=sub->entity[i].scvrt;		 temp_el[k].scvrt=map[n2];
	n3=sub->entity[i].thvrt;		 temp_el[k].thvrt=map[n3];
	srf=sub->entity[i].ms_wert; temp_el[k].ms_wert=srf;
	k++;
	}
sub->e_grs=k;

for(i=0;i<sub->n_grs;i++)
	{sub->knot[i].absi=temp_nd[i].absi;
	sub->knot[i].ordo=temp_nd[i].ordo;
	sub->knot[i].cote=temp_nd[i].cote;
	}

for(i=0;i<sub->e_grs;i++)
	{sub->entity[i].frvrt=temp_el[i].frvrt;
	sub->entity[i].scvrt=temp_el[i].scvrt;
	sub->entity[i].thvrt=temp_el[i].thvrt;
	sub->entity[i].ms_wert=temp_el[i].ms_wert;
	}
free(temp_el);
free(temp_nd);
free(map);
}



double qosc_anis_vord(manif_tl sub)
{int e,nel,n1,n2,n3;
double ar,loc;
nel=sub.e_grs;
ar=1.0;
for(e=0;e<nel;e++)
	{n1=sub.entity[e].frvrt;
	n2=sub.entity[e].scvrt;
	n3=sub.entity[e].thvrt;
	loc=potb_anis_dopg(sub.knot[n1],sub.knot[n2],sub.knot[n3]);
	if(loc<ar)
		ar=loc;
	}
return ar;
}



void tugw_extr_vuzf(manif_tl msh,manif_tl *sub,int *elm,int nb_sub,
int *map_node,int *map_elem,int max_node,int max_elem)
{int i,j,k,*map,*inv,nnd,s,n[4],ts,dummy,m;
nnd=msh.n_grs;
map=(int *)malloc(3*nb_sub*sizeof(int));
inv=(int *)malloc(nnd*sizeof(int));

k=0;
for(i=0;i<nb_sub;i++)
	{s=elm[i];
	n[1]=msh.entity[s].frvrt;
	n[2]=msh.entity[s].scvrt;
	n[3]=msh.entity[s].thvrt;
	for(j=1;j<=3;j++)
		{ts=gonl_arra_govj(map,k,n[j],&dummy);
		if(ts==0)
			{map[k]=n[j];
			inv[n[j]]=k;
			k++;
			}
		}
	}

sub->n_grs=k;
for(i=0;i<sub->n_grs;i++)
	{s=map[i];
	if(i>=max_node)
		{fprintf(tmpout,"max_node is reached\n");
		exit(0);
		}
	getf_find_rogc_todj(msh.knot[s],&sub->knot[i]);
	map_node[i]=s;
	}

sub->e_grs=nb_sub;
for(i=0;i<nb_sub;i++)
	{s=elm[i];
	if(i>=max_elem)
		{fprintf(tmpout,"max_elem[%d] is reached\n",max_elem);
		exit(0);
		}
	m=msh.entity[s].frvrt;	sub->entity[i].frvrt=inv[m];
	m=msh.entity[s].scvrt;	sub->entity[i].scvrt=inv[m];
	m=msh.entity[s].thvrt;	sub->entity[i].thvrt=inv[m];
	sub->entity[i].ms_wert=msh.entity[s].ms_wert;
	map_elem[i]=s;
	}
free(map);
free(inv);
}



double nich_pros_fiph(manif_tl *sub,int e,point P)
{int m1,m2;
double ar;
m1=sub->kt[e].frvrt;
m2=sub->kt[e].scvrt;
folp_fuse_ruqk(sub,m1,m2,e,P);
ar=qosc_anis_vord(*sub);
return ar;
}



int dils_find_fisv(manif_tl sub,int n1,int n2,int *map_node)
{int e=-1,i,ned,m1,m2,N1,N2;
ned=sub.k_grs;
for(i=0;i<ned;i++)
	{m1=sub.kt[i].frvrt;
	m2=sub.kt[i].scvrt;
	N1=map_node[m1];
	N2=map_node[m2];
	if(((N1==n1)&&(N2==n2))||((N1==n2)&&(N2==n1)))
		{e=i;
		break;
		}
	}
if(e==-1)
	{fprintf(tmpout,"1-Unable to find local edge\n");
	exit(0);
	}
return e;
}



int gubs_feas_dujz(manif_tl msh,int e,double anis,point P)
{int *elm,i,j,k,max_elm=200,nd[2],m,val,ts,res,e_loc,n1,n2;
int *temp,dummy,nb_sub,*map_node,*map_elem,max_ed_loc;
double ar;
manif_tl sub;
elm=(int *)malloc(max_elm*sizeof(int));
nd[0]=msh.kt[e].frvrt;
nd[1]=msh.kt[e].scvrt;

k=0;
for(i=0;i<2;i++)
	{m=nd[i];
	val=msh.increl[m].val;
	temp=(int *)malloc(2*val*sizeof(int));
	cepj_inci_jeqt(msh,m,temp);
	for(j=0;j<val;j++)
		{ts=gonl_arra_govj(elm,k,temp[j],&dummy);
		if(ts==0)
			{if(k>=max_elm)
				{fprintf(tmpout,"max_elm is reached\n");
				exit(0);
				}
			elm[k]=temp[j];
			k++;
			}
		}
	free(temp);
	}
nb_sub=k;

sub.entity=(telolf *)malloc(nb_sub*sizeof(telolf));
sub.knot   =(point *)malloc(3*nb_sub*sizeof(point));
map_node   =(int *)malloc(3*nb_sub*sizeof(int));
map_elem   =(int *)malloc(nb_sub*sizeof(int));
tugw_extr_vuzf(msh,&sub,elm,nb_sub,map_node,map_elem,3*nb_sub,nb_sub);

max_ed_loc=3*sub.e_grs;
sub.kt=(kt *)malloc(max_ed_loc*sizeof(kt));
cogv_fill_zicd(&sub,max_ed_loc);
n1=msh.kt[e].frvrt;
n2=msh.kt[e].scvrt;
e_loc=dils_find_fisv(sub,n1,n2,map_node);
ar=nich_pros_fiph(&sub,e_loc,P);
if(ar>anis)
	res=SUCCESS;
else
	res=FAILURE;
free(elm);
free(sub.kt);
free(sub.entity);
free(sub.knot);
free(map_node);
free(map_elem);
return res;
}


void vulp_find_rahk_vork(fund_quad Q_in,
fund_quad *Q_out)
{int i,j;
for(i=0;i<4;i++)
for(j=0;j<4;j++)
	Q_out->A[i][j]=Q_in.A[i][j];
for(i=0;i<3;i++)
	Q_out->b[i]=Q_in.b[i];
Q_out->c=Q_in.c;
}


int cefg_find_fovw(manif_tl msh,int n1,int n2)
{int val,i,e,n_a,n_b,res=-1;
val=msh.increl[n1].val;
for(i=0;i<val;i++)
	{e=msh.increl[n1].inc[i];
	n_a=msh.kt[e].frvrt;
	n_b=msh.kt[e].scvrt;
	if((n_a==n2)||(n_b==n2))
		{res=e;
		break;
		}
	}
return res;
}


int goqv_inci_liqf(manif_tl msh,int e,int *inv_node,
int n_inv,int *map_nd)
{int nd[2],i,j,k,val,E,m1,m2,md,*inc,val0,val1;
int ts,dummy,res;
nd[0]=msh.kt[e].frvrt;
nd[1]=msh.kt[e].scvrt;
val0=msh.increl[nd[0]].val;
val1=msh.increl[nd[1]].val;
inc=(int *)malloc((val0+val1)*sizeof(int));
k=0;
for(i=0;i<2;i++)
	{val=msh.increl[nd[i]].val;
	for(j=0;j<val;j++)
		{E=msh.increl[nd[i]].inc[j];
		m1=msh.kt[E].frvrt;
		m2=msh.kt[E].scvrt;
		if(m1==nd[i])	md=m2;
		if(m2==nd[i])	md=m1;
		ts=gonl_arra_govj(inc,k,md,&dummy);
		if(ts==0)
			{inc[k]=md;
			k++;
			}
		}
	}

res=0;
for(i=0;i<k;i++)
	{ts=gonl_arra_govj(inv_node,n_inv,inc[i],&dummy);
	if(ts==1)
		{res=1;
		break;
		}
	}
free(inc);
return res;
}


void kijn_veri_jorf(double *err,fund_quad *Q,manif_tl msh)
{int e,ned,n1,n2;
double er,diff;
point dummy;
ned=msh.k_grs;
for(e=0;e<ned;e++)
	{n1=msh.kt[e].frvrt;
	n2=msh.kt[e].scvrt;
	er=selr_opti_zoqp(msh,n1,n2,Q[n1],Q[n2],&dummy);
	diff=fabs(er-err[e]);
	if(diff>1.0e-8)
		{fprintf(tmpout,"Dissimilar errors\n");
		fprintf(tmpout,"err1=%e   err2=%e\n",er,err[e]);
		exit(0);
		}
	}
fprintf(tmpout,"GOOD QUADRICS EDGE ERRORS\n");
}


void geqw_disp_kudv(fund_quad Q)
{int i;
fprintf(tmpout,"A:\n");
for(i=0;i<3;i++)
	fprintf(tmpout,"%f   %f   %f\n",Q.A[i][0],Q.A[i][1],Q.A[i][2]);
fprintf(tmpout,"b=[%f,%f,%f]  c=%f\n",Q.b[0],Q.b[1],Q.b[2],Q.c);
}


double werg_dist_fenv(fund_quad Q1,
fund_quad Q2)
{int i,j;
double res;
res=0.0;
for(i=0;i<3;i++)
	{for(j=0;j<3;j++)
		res=res+fabs(Q1.A[i][j]-Q2.A[i][j]);
	res=res+fabs(Q1.b[i]-Q2.b[i]);
	}
res=res+fabs(Q1.c-Q2.c);
return res;
}


void guql_veri_qacj(fund_quad *Q_list,manif_tl msh)
{int z,nnd;
double dis;
fund_quad Q;
fprintf(tmpout,"Verify nodal errors\n");
nnd=msh.n_grs;
for(z=0;z<nnd;z++)
	{qemt_dire_lajk(msh,z,&Q);
	dis=werg_dist_fenv(Q,Q_list[z]);
	if(dis>0.01)
		{geqw_disp_kudv(Q);
		geqw_disp_kudv(Q_list[z]);
		exit(0);
		}
	}
fprintf(tmpout,"GOOD QUADRICS NODAL ERRORS\n");
}


void jofk_veri_corg(manif_tl msh)
{int nel,i,n1,n2,n3;
double ar1,ar2,diff,eps=1.0e-4;
nel=msh.e_grs;
for(i=0;i<nel;i++)
	{n1 =msh.entity[i].frvrt;
	n2  =msh.entity[i].scvrt;
	n3  =msh.entity[i].thvrt;
	ar1 =valp_area_qelk(msh.knot[n1],msh.knot[n2],msh.knot[n3]);
	ar2 =msh.entity[i].ms_wert;
	diff=fabs(ar1-ar2);
	if(diff>eps)
		{fprintf(tmpout,"Verify measure entry\n");
		exit(0);
		}
	}
fprintf(tmpout,"GOOD MEASURE ENTRY\n");
}



void ruqn_upda_ducz(manif_tl *msh,int *old_n1,int *old_n2,
int *inv_node,int n_inv,int *map_nd,int nnd_old,int ned_old,
double *err,fund_quad *Q)
{int i,ts1,ts2,im1,im2,e_new,z;
int n1,n2,ts,dummy,ts_inc;
double *err_temp;
fund_quad *Q_temp;
point P;

if(msh->k_grs>ned_old)
	{fprintf(tmpout,"WARNING: increasing kt count\n");
	exit(0);
	}
err_temp=(double *)malloc(msh->k_grs*sizeof(double));
for(i=0;i<msh->k_grs;i++)
	err_temp[i]=CHANGE_ID;
for(i=0;i<ned_old;i++)
	{n1=old_n1[i];
	n2=old_n2[i]; 
	ts1=gonl_arra_govj(inv_node,n_inv,n1,&dummy);
	if(ts1==0)
		{ts2=gonl_arra_govj(inv_node,n_inv,n2,&dummy);
		if(ts2==0)
			{im1=map_nd[n1];
			im2=map_nd[n2]; 
			e_new=cefg_find_fovw(*msh,im1,im2);
			if(e_new==-1)
				{fprintf(tmpout,"Unable to find e_new\n");
				fprintf(tmpout,"[n1,n2]=[%d,%d]\n",n1,n2);
				fprintf(tmpout,"[im1,im2]=[%d,%d]\n",im1,im2);
				exit(0);
				}
			ts_inc=goqv_inci_liqf(*msh,e_new,inv_node,n_inv,map_nd);
			if(ts_inc==0)
				err_temp[e_new]=err[i];
			}
		}
	}
for(i=0;i<msh->k_grs;i++)
	err[i]=err_temp[i];
free(err_temp);

Q_temp=(fund_quad *)malloc(nnd_old*sizeof(fund_quad));
for(i=0;i<nnd_old;i++)
	{z=map_nd[i];
	ts=gonl_arra_govj(inv_node,n_inv,i,&dummy);
	if(ts==0)
		vulp_find_rahk_vork(Q[i],&Q_temp[z]);
	else
		{if(z!=-1)
			qemt_dire_lajk(*msh,z,&Q_temp[z]);
		}
	}



for(i=0;i<msh->n_grs;i++)
	vulp_find_rahk_vork(Q_temp[i],&Q[i]);
free(Q_temp);

for(i=0;i<msh->k_grs;i++)
if(fabs(err[i]-CHANGE_ID)<1.0e-7)
	{n1=msh->kt[i].frvrt;
	n2=msh->kt[i].scvrt;
	err[i]=selr_opti_zoqp(*msh,n1,n2,Q[n1],Q[n2],&P);
	}

}

 

int qulw_vois_nuhk(manif_tl msh,int n1,int n2,int *nd_vs)
{int nd[2],i,j,nb,z,nd_a,nd_b,val,E;
int vs,ts,dummy;
nd[0]=n1;	nd[1]=n2;
nb=0;
for(i=0;i<2;i++)
	{z=nd[i];
	val=msh.increl[z].val;
	for(j=0;j<val;j++)
		{E=msh.increl[z].inc[j];
		nd_a=msh.kt[E].frvrt;
		nd_b=msh.kt[E].scvrt;
		if(nd_a==z)	vs=nd_b;
		if(nd_b==z)	vs=nd_a;
		ts=gonl_arra_govj(nd_vs,nb,vs,&dummy);
		if(ts==0)
			{nd_vs[nb]=vs;
			nb++;
			}
		}
	}
return nb;
}


