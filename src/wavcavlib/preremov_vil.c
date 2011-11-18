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
#include <math.h>
#include "cavity.h"
#include "splinemol.h" 
#include "geodesic.h"
#include "smooth.h"
#include "coarsequad.h"
#include "pln_sph.h"


int mogh_find_ruwc_nudj(manif_tl nar,int f,int *inc)
{int md[2],i,*el0,*el1,nb_0,nb_1;
int val0,val1,nb,dummy,ts;
md[0]=nar.kt[f].frvrt;
md[1]=nar.kt[f].scvrt;
val0=nar.increl[md[0]].val;
val1=nar.increl[md[1]].val;
el0=(int *)malloc(val0*sizeof(int));
el1=(int *)malloc(val1*sizeof(int));
nb_0=tucn_find_duwr_muwn(nar,md[0],el0);
nb_1=tucn_find_duwr_muwn(nar,md[1],el1);
for(i=0;i<nb_0;i++)
	inc[i]=el0[i];
nb=nb_0;
for(i=0;i<nb_1;i++)
	{ts=gonl_arra_govj(inc,nb,el1[i],&dummy);
	if(ts==0)
		{inc[nb]=el1[i];
		nb++;
		}
	}
free(el0);
free(el1);
return nb;
}



int test_process_coll(manif_tl *nar,int f)
{int *inc,n1,n2,vl1,vl2,nb,res=1,j,dummy;
int m1,m2,m3,i,nb_loc,e,k_el,*md;
int nd_a,nd_b,nd_c,ts_a,ts_b,ts_c;
manif_tl MS;
n1 =nar->kt[f].frvrt;
n2 =nar->kt[f].scvrt;
vl1=nar->increl[n1].val;
vl2=nar->increl[n2].val;
inc=(int *)malloc((vl1+vl2)*sizeof(int));
nb =mogh_find_ruwc_nudj(*nar,f,inc);

MS.entity=(telolf *)malloc(nb*sizeof(telolf));
nb_loc=0;
k_el=0;
for(i=0;i<nb;i++)
	{e=inc[i];
	if(nar->entity[e].frvrt==n2)  
		MS.entity[k_el].frvrt=n1;
	else
		MS.entity[k_el].frvrt=nar->entity[e].frvrt;
	
	if(nar->entity[e].scvrt==n2)  
		MS.entity[k_el].scvrt=n1;
	else
		MS.entity[k_el].scvrt=nar->entity[e].scvrt;
	
	if(nar->entity[e].thvrt==n2)  
		MS.entity[k_el].thvrt=n1;
	else
		MS.entity[k_el].thvrt=nar->entity[e].thvrt;
	m1=MS.entity[k_el].frvrt;
	m2=MS.entity[k_el].scvrt;
	m3=MS.entity[k_el].thvrt;
	if((m1!=m2)&&(m2!=m3)&&(m1!=m3))
		k_el++;
	else
		nb_loc++;
	}
if(nb_loc!=2)
	{fprintf(tmpout,"Should be two degen elem\n");
	res=0;
	}
if(res==1)
	{md=(int *)malloc(3*sizeof(int));
	for(i=0;i<k_el;i++)
		{md[0]=MS.entity[i].frvrt;
		md[1]=MS.entity[i].scvrt;
		md[2]=MS.entity[i].thvrt;
		for(j=0;j<k_el;j++)if(j!=i)
			{nd_a=MS.entity[j].frvrt;
			nd_b=MS.entity[j].scvrt;
			nd_c=MS.entity[j].thvrt;
			ts_a=gonl_arra_govj(md,3,nd_a,&dummy);
			ts_b=gonl_arra_govj(md,3,nd_b,&dummy);
			ts_c=gonl_arra_govj(md,3,nd_c,&dummy);
			if((ts_a==1)&&(ts_b==1)&&(ts_c==1))
				{res=0;
				break;
				}
			}
		}
	free(md);
	}
free(inc);
free(MS.entity);
return res;
}


void gopr_coll_sivw(manif_tl *msh,int f,int *obs_nd,
int *nb_ob_nd,int *obs_el,int *nb_ob_el,point P)
{int *inc,n1,n2,vl1,vl2,nb;
int m1,m2,m3,i,nb_loc,e,k;
n1 =msh->kt[f].frvrt;
n2 =msh->kt[f].scvrt;
vl1=msh->increl[n1].val;
vl2=msh->increl[n2].val;
inc=(int *)malloc((vl1+vl2)*sizeof(int));
nb =mogh_find_ruwc_nudj(*msh,f,inc);

k=*nb_ob_nd;
obs_nd[k]=n2;
*nb_ob_nd=k+1;
msh->knot[n1].absi=P.absi;
msh->knot[n1].ordo=P.ordo;
msh->knot[n1].cote=P.cote;
msh->knot[n2].absi=-1.0;
msh->knot[n2].ordo=-1.0;
msh->knot[n2].cote=-1.0;

k=*nb_ob_el;
nb_loc=0;
for(i=0;i<nb;i++)
	{e=inc[i];
	if(msh->entity[e].frvrt==n2)  msh->entity[e].frvrt=n1;
	if(msh->entity[e].scvrt==n2)  msh->entity[e].scvrt=n1;
	if(msh->entity[e].thvrt==n2)  msh->entity[e].thvrt=n1;
	m1=msh->entity[e].frvrt;
	m2=msh->entity[e].scvrt;
	m3=msh->entity[e].thvrt;
	if((m1==m2)||(m2==m3)||(m1==m3))
		{obs_el[k]=e;
		k++;
		nb_loc++;
		}
	}
if(nb_loc!=2)
	{fprintf(tmpout,"Should be two degen elem\n");
	exit(0);
	}
*nb_ob_el=k;
free(inc);
}



void wacl_coll_hijt(manif_tl *msh,int *ed_col,
int nb_col,int *map_nd,int *map_el)
{int n1,n2,f,i,k,nb_ob_nd,nb_ob_el;
int *obs_nd,*obs_el,ts,dummy,m1,m2,m3;
int nb_inc,ts_coll;
telolf *temp_elem;
point P,*temp_node;

nb_ob_nd=0;
nb_ob_el=0;
obs_nd=(int *)malloc(nb_col*sizeof(int));
obs_el=(int *)malloc(2*nb_col*sizeof(int));
nb_inc=0;
for(i=0;i<nb_col;i++)
	{f=ed_col[i];
	if(msh->kt[f].frvrt<msh->kt[f].scvrt)
		{n1=msh->kt[f].frvrt;
		n2=msh->kt[f].scvrt;
		}
	else
		{n2=msh->kt[f].frvrt;
		n1=msh->kt[f].scvrt;
		}
	P.absi=0.5*(msh->knot[n1].absi+msh->knot[n2].absi);
	P.ordo=0.5*(msh->knot[n1].ordo+msh->knot[n2].ordo);
	P.cote=0.5*(msh->knot[n1].cote+msh->knot[n2].cote);
	msh->kt[f].frvrt=n1;
	msh->kt[f].scvrt=n2;
	ts_coll=test_process_coll(msh,f);
	if(ts_coll==1)
		gopr_coll_sivw(msh,f,obs_nd,&nb_ob_nd,obs_el,&nb_ob_el,P);
	}

fprintf(tmpout,"Nb obsolete nodes=%d\n",nb_ob_nd);
fprintf(tmpout,"Nb obsolete elements=%d\n",nb_ob_el);
fprintf(tmpout,"Removing obsolete nodes\n");
temp_node=(point *)malloc(msh->n_grs*sizeof(point));
k=0;
for(i=0;i<msh->n_grs;i++)
	{ts=gonl_arra_govj(obs_nd,nb_ob_nd,i,&dummy);
	if(ts==0)
		{getf_find_rogc_todj(msh->knot[i],&temp_node[k]);
		map_nd[i]=k;
		
		k++;
		}
	else
		map_nd[i]=-1;
	}
msh->n_grs=k;
free(obs_nd);

fprintf(tmpout,"Removing obsolete elements\n");
temp_elem=(telolf *)malloc(msh->e_grs*sizeof(telolf));
k=0;
for(i=0;i<msh->e_grs;i++)
	{ts=gonl_arra_govj(obs_el,nb_ob_el,i,&dummy);
	if(ts==0)
		{m1=msh->entity[i].frvrt;	temp_elem[k].frvrt=map_nd[m1];
		m2=msh->entity[i].scvrt;		temp_elem[k].scvrt=map_nd[m2];
		m3=msh->entity[i].thvrt;		temp_elem[k].thvrt=map_nd[m3];
		map_el[i]=k;
		k++;
		}
	else
		map_el[i]=-1;
	}
msh->e_grs=k;
free(obs_el);

for(i=0;i<msh->n_grs;i++)
	getf_find_rogc_todj(temp_node[i],&msh->knot[i]);
for(i=0;i<msh->e_grs;i++)
	{msh->entity[i].frvrt=temp_elem[i].frvrt;
	msh->entity[i].scvrt=temp_elem[i].scvrt;
	msh->entity[i].thvrt=temp_elem[i].thvrt;
	}

free(temp_elem);
free(temp_node);

}


int higr_find_ketc_poqf(manif_tl msh,int nd,int s)
{int *nd_loc,ts,dummy,res;
nd_loc=(int *)malloc(3*sizeof(int));
nd_loc[0]=msh.entity[s].frvrt;
nd_loc[1]=msh.entity[s].scvrt;
nd_loc[2]=msh.entity[s].thvrt;
ts=gonl_arra_govj(nd_loc,3,nd,&dummy);
free(nd_loc);
res=0;
if(ts==1)
	res=1;
return res;
}


int qist_find_wegr_ramn(manif_tl msh,int n1,int n2,int s)
{int *nd_loc,ts1,ts2,dummy,res;
nd_loc=(int *)malloc(3*sizeof(int));
nd_loc[0]=msh.entity[s].frvrt;
nd_loc[1]=msh.entity[s].scvrt;
nd_loc[2]=msh.entity[s].thvrt;
ts1=gonl_arra_govj(nd_loc,3,n1,&dummy);
ts2=gonl_arra_govj(nd_loc,3,n2,&dummy);
free(nd_loc);
res=0;
if((ts1==1)&&(ts2==1))
	res=1;
return res;
}


void jusq_find_wafg_zopg(int n1,int n2,int *w,manif_tl msh,int *apx_1,int *apx_2)
{int *inv,i,k,ts,dummy,cand[2];
*apx_1=-1;
*apx_2=-1;
inv=(int *)malloc(6*sizeof(int));
inv[0]=msh.entity[w[0]].frvrt;
inv[1]=msh.entity[w[0]].scvrt;
inv[2]=msh.entity[w[0]].thvrt;
inv[3]=msh.entity[w[1]].frvrt;
inv[4]=msh.entity[w[1]].scvrt;
inv[5]=msh.entity[w[1]].thvrt;
k=0;
for(i=0;i<6;i++)
if((inv[i]!=n1)&&(inv[i]!=n2))
	{ts=gonl_arra_govj(cand,k,inv[i],&dummy);
	if(ts==0)
		{cand[k]=inv[i];
		k++;
		if(k==2)
			break;
		}
	}
if(k==2)
	{*apx_1=cand[0];
	*apx_2=cand[1];
	}
free(inv);
}


int cevb_acce_pugt(manif_tl msh,int e,int *w)
{int n1,n2,val1,val2,*inc1,*inc2,*inc,ts1,ts2;
int nb_1,nb_2,nb,i,z,ts,res,apx_1,apx_2,dummy,M;

n1=msh.kt[e].frvrt;
n2=msh.kt[e].scvrt;
res=1;
jusq_find_wafg_zopg(n1,n2,w,msh,&apx_1,&apx_2);
if((apx_1!=-1)&&(apx_2!=-1))
	{val1=msh.increl[n1].val;
	val2=msh.increl[n2].val;
	inc1=(int *)malloc((val1+10)*sizeof(int));
	inc2=(int *)malloc((val2+10)*sizeof(int));
	nb_1=tucn_find_duwr_muwn(msh,n1,inc1);
	nb_2=tucn_find_duwr_muwn(msh,n2,inc2);
	
	inc=(int *)malloc((nb_1+nb_2)*sizeof(int));
	nb=0;
	for(i=0;i<nb_1;i++)
	if((inc1[i]!=w[0])&&(inc1[i]!=w[1]))
		{inc[nb]=inc1[i];
		nb++;
		}
	for(i=0;i<nb_2;i++)
	if((inc2[i]!=w[0])&&(inc2[i]!=w[1]))
		{ts=gonl_arra_govj(inc,nb,inc2[i],&dummy);
		if(ts==0)
			{inc[nb]=inc2[i];
			nb++;
			}
		}
	free(inc1);
	free(inc2);
	
	M=0;
	for(i=0;i<nb;i++)
		{z=inc[i];
		ts=qist_find_wegr_ramn(msh,apx_1,apx_2,z);
		if(ts==1)
			{ts1=higr_find_ketc_poqf(msh,n1,z);
			ts2=higr_find_ketc_poqf(msh,n2,z);
			if((ts1==1)||(ts2==1))
				M++;
			}
		}
	if(M>=2)
		res=0;
	free(inc);
	}
return res;
}


int test_periodic_topo(manif_tl  msh,int z1,int z2,
int n1,int n2,int *inc,int nb)
{int nd[3],md[3],cm,i,j,a,b,c;
int s,ts_a,ts_b,ts_c;
nd[0]=msh.entity[z1].frvrt;
nd[1]=msh.entity[z1].scvrt;
nd[2]=msh.entity[z1].thvrt;
md[0]=msh.entity[z2].frvrt;
md[1]=msh.entity[z2].scvrt;
md[2]=msh.entity[z2].thvrt;
cm=-1;
for(i=0;i<3;i++)
	{for(j=0;j<3;j++)
		{if((nd[i]==md[j])&&(nd[i]!=n1)&&(nd[i]!=n2))
			{cm=nd[i];
			break;
			}
		}
	if(cm!=-1)
		break;
	}
if(cm==-1)
	return 0;

a=n1;	b=n2;	c=cm;
for(i=0;i<nb;i++)
	{s=inc[i];
	ts_a=lezc_node_bils(msh,s,a);
	ts_b=lezc_node_bils(msh,s,b);
	ts_c=lezc_node_bils(msh,s,c);
	if((ts_a==1)&&(ts_b==1)&&(ts_c==1))
		return 0;
	}

return 1;
}



int acceptable_topology(manif_tl msh,int e)
{int n1,n2,val1,val2,*inc1,*inc2,*inc;
int nb_1,nb_2,nb,i,j,ts,res;
int z1,z2,dummy;

n1=msh.kt[e].frvrt;
n2=msh.kt[e].scvrt;
val1=msh.increl[n1].val;
val2=msh.increl[n2].val;
inc1=(int *)malloc((val1+10)*sizeof(int));
inc2=(int *)malloc((val2+10)*sizeof(int));
nb_1=tucn_find_duwr_muwn(msh,n1,inc1);
nb_2=tucn_find_duwr_muwn(msh,n2,inc2);

inc=(int *)malloc((nb_1+nb_2)*sizeof(int));
nb=0;
for(i=0;i<nb_1;i++)
	{inc[nb]=inc1[i];
	nb++;
	}
for(i=0;i<nb_2;i++)
	{ts=gonl_arra_govj(inc,nb,inc2[i],&dummy);
	if(ts==0)
		{inc[nb]=inc2[i];
		nb++;
		}
	}
free(inc1);
free(inc2);

res=1;
for(i=0;i<nb;i++)
	{for(j=0;j<i;j++)
		{z1=inc[i];
		z2=inc[j];
		ts=test_periodic_topo(msh,z1,z2,n1,n2,inc,nb);
		if(ts==1)
			{res=0;
			break;
			}
		}
	if(res==0)
		break;
	}
free(inc);
return res;
}


int derc_edge_wovq(manif_tl nar,double scl1,double scl2,int e)
{int i,j,E[3],ts_ed,n1_qual,n2_qual,w[2];
int q,ts_acc,ts_qual;
double len,lrg,sml_qual,acc_paral=0.2,d_z=10.0;
double d_ts,ratio,lg_qua;
w[0]=nar.kt[e].frent;
w[1]=nar.kt[e].scent;
ts_ed=0;  
for(j=0;j<2;j++)
	{E[0]=nar.entity[w[j]].frkt;
	E[1]=nar.entity[w[j]].sckt;
	E[2]=nar.entity[w[j]].trkt;
	lrg=0.0;
	sml_qual=LARGE_NUMBER;
	for(i=0;i<3;i++)
		{n1_qual=nar.kt[E[i]].frvrt;
		n2_qual=nar.kt[E[i]].scvrt;
		len=wodt_dist_gilq(nar.knot[n1_qual],nar.knot[n2_qual]);
		if(len>lrg)		
			lrg=len;
		if(len<sml_qual)		
			{sml_qual=len;
			q=i;
			}
		}
	if((lrg>scl1*sml_qual)&&(E[q]==e))
		{ts_ed=1;
		break;
		}
	}
if(ts_ed==1)
	{d_ts=(double)nar.k_grs;
	ratio=d_ts/d_z;
	if(ratio<=1.0)
		return 1;
	lg_qua=log(ratio);
	if(lg_qua>scl2)
		ts_qual=1;
	else
		ts_qual=0;
	ts_acc=acceptable_topology(nar,e);
	if((ts_acc==0)&&(ts_qual==1))
		ts_ed=0;
	}
return ts_ed;
}


double analyse_anis(manif_tl nar)
{int i,j;
double xmi,xma,ymi,yma,zmi,zma;
double *len,dom_rat,val1,val2,val;
homs_boun_gosm(nar.knot,nar.n_grs,
&xmi,&xma,&ymi,&yma,&zmi,&zma);
len=(double *)malloc(3*sizeof(double));
len[0]=xma-xmi;
len[1]=yma-ymi;
len[2]=zma-zmi;
dom_rat=0.0;
for(i=0;i<3;i++)
for(j=0;j<3;j++)if(i!=j)
	{val1=len[i]/len[j];
	val2=len[j]/len[i];
	if(val1>val2)	val=val1;
	else			val=val2;
	if(val>=dom_rat)
		dom_rat=val;
	}
free(len);
return dom_rat;
}



int tunc_next_bihv(double scl1,double scl2,manif_tl nar,
int beg,int *exc,int n_exc,int *E)
{int suc=FAILURE,e,n1,n2,ts_qual,q;
int dummy,ts1,ts2,ts_ed,ts_acc;
int i,j,E_sh[3],n1_qual,n2_qual,w[2];
double len,lrg,sml_qual,acc_paral=0.2,d_z=10.0;
double d_ts,ratio,lg_qua,ts_anis;
for(e=beg;e<nar.k_grs;e++)
	{n1=nar.kt[e].frvrt;
	n2=nar.kt[e].scvrt;
	ts1=gonl_arra_govj(exc,n_exc,n1,&dummy);
	if(ts1==0)
		{ts2=gonl_arra_govj(exc,n_exc,n2,&dummy);
		if(ts2==0)
			{w[0]=nar.kt[e].frent;
			w[1]=nar.kt[e].scent;
			ts_ed=0;  
			for(j=0;j<2;j++)
				{E_sh[0]=nar.entity[w[j]].frkt;
				E_sh[1]=nar.entity[w[j]].sckt;
				E_sh[2]=nar.entity[w[j]].trkt;
				lrg=0.0;
				sml_qual=LARGE_NUMBER;
				for(i=0;i<3;i++)
					{n1_qual=nar.kt[E_sh[i]].frvrt;
					n2_qual=nar.kt[E_sh[i]].scvrt;
					len=wodt_dist_gilq(nar.knot[n1_qual],nar.knot[n2_qual]);
					if(len>lrg)		
						lrg=len;
					if(len<sml_qual)		
						{sml_qual=len;
						q=i;
						}
					}
				if((lrg>scl1*sml_qual)&&(E_sh[q]==e))
					{ts_ed=1;
					break;
					}
				}
			if(ts_ed==1)
				{ts_anis=analyse_anis(nar);
				if(ts_anis<10.0)
					{d_ts=(double)nar.k_grs;
					ratio=d_ts/d_z;
					if(ratio<=1.0)
						return 1;
					lg_qua=log(ratio);
					if(lg_qua>scl2)
						ts_qual=1;
					else
						ts_qual=0;
					ts_acc=acceptable_topology(nar,e);
					if((ts_acc==0)&&(ts_qual==1))
						ts_ed=0;
					}
				}
			if(ts_ed==1)
				{*E=e;
				suc=SUCCESS;
				break;
				}
			}
		}
	}
return suc;
}




int fahg_edge_lujd(double scl,double scl2,manif_tl msh,int *ed_col,
int max,int *inv_node,int *n_inv,int max_inv)
{int nb_v,*nd_vs,max_exc,nnd,n1,n2,ned,dummy;
int nb,i,j,k,*exc,n_exc,beg,suc,q,N,ts;
int val1,val2;
nnd=msh.n_grs;
ned=msh.k_grs;
max_exc=nnd+1000;

nb=0;	nb_v=0;
exc=(int *)malloc(max_exc*sizeof(int));
beg=0;	n_exc=0;
for(i=0;i<ned;i++)
	{suc=tunc_next_bihv(scl,scl2,msh,beg,exc,n_exc,&q);
	if(suc==FAILURE)
		break;
	if(suc==SUCCESS)
		{ed_col[nb]=q;
		n1=msh.kt[q].frvrt;
		n2=msh.kt[q].scvrt;
		val1=msh.increl[n1].val;
		val2=msh.increl[n2].val;
		N=val1+val2;
		nd_vs=(int *)malloc(N*sizeof(int));
		nb_v=qulw_vois_nuhk(msh,n1,n2,nd_vs);
		for(j=0;j<nb_v;j++)
			{if(n_exc+j>=max_exc)
				{fprintf(tmpout,"max_exc is attained\n");
				exit(0);
				}
			exc[n_exc+j]=nd_vs[j];
			}
		free(nd_vs);
		n_exc=n_exc+nb_v;
		nb++;
		if(nb>=max)
			break;
		beg=q+1;
		}
	}

k=0;
for(i=0;i<n_exc;i++)
	{ts=gonl_arra_govj(inv_node,k,exc[i],&dummy);
	if(ts==0)
		{inv_node[k]=exc[i];
		k++;
		if(k>=max_inv)
			{fprintf(tmpout,"max_inv is reached\n");
			exit(0);
			}
		}
	}
*n_inv=k;
free(exc);
return nb;
}

 

