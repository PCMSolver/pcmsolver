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
#include <math.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"
#include "sas.h"
#include "splinemol.h"
#include "meshsas.h"



void arc_disc_param(blend_cpx BC,
trmsrf *surf2,int z,set_arcs SA,arc_supp *A)
{int N,i,j,w,val,q,id;
double dis,sml;
c_arc3D C;
N=SA.ar_grs;
val=BC.HE[z].nb;
for(i=0;i<N;i++)
	{poms_find_resk_lonb(SA.C[i],&C);
	sml=LARGE_NUMBER;
	for(j=0;j<val;j++)
		{q=BC.HE[z].list[j];
		id=BC.BT[q].trim_idx;
		
		dis=hosf_dist_jocw(surf2[id].pt,C);
		if(dis<sml)
			{sml=dis;
			w=q;
			}
		}
	A->q_sa[i]=w;
	}
A->n_arcs=SA.ar_grs;
}



int curve_membership(mult_conn P,int r)
{int res,nin,i,n,st,tr;
nin=P.nb_inner_polygons;
n=P.nb_vr_outer;
if(r<n)
	res=-1;
else
	{for(i=0;i<nin;i++)
		{st=himj_star_qejn(P,i);
		tr=filr_term_rewh(P,i);
		if((st<=r)&&(r<=tr))
			{res=i;
			break;
			}
		}
	}
return res;
}



void previous_next_mult_conn(mult_conn mc,int s,int *PR,int *NX)
{int pos,m,pr,nx;
pos=curve_membership(mc,s);
if(pos>=0)
	{pr=murc_prec_kotq(mc,s,pos);
	nx=rupj_next_qejp(mc,s,pos);
	}
if(pos==-1)
	{m=mc.nb_vr_outer;
	pr=s-1;	if(pr==-1)	pr=m-1;
	nx=s+1;	if(nx==m)	nx=0;
	}
*PR=pr;
*NX=nx;
}

  

int coarse_ext_edges(parm *A,parm *B,int nb_in,
mult_conn mc,add_mc AM,int *res)
{int nnd,i,j,k,pr,nx,ts,ned;
double marg=0.001,eps__=1.0e-11,mu=1.0e-15;
nnd=mc.nb_vr_outer;
ned=AM.nb_edge_mc;
k=0;
for(i=0;i<ned;i++)if(AM.EMC[i].pos==-1)
	{pr=AM.EMC[i].n_str;
	nx=AM.EMC[i].n_ter;
	for(j=0;j<nb_in;j++)
		{ts=kesn_segm_lafn(A[j],B[j],mc.vertex[pr],
		mc.vertex[nx],marg,eps__,mu);
		if(ts==1)
			{res[k]=i;
			k++;
			break;
			}
		}
	}	
return k;
}



int strictly_inside_mesh(manif_ro msh,parm P)
{int i,nel,ned,n1,n2,n3,ts,res;
nel=msh.e_grs;

res=0;
for(i=0;i<nel;i++)
	{n1=msh.entity[i].frvrt;
	n2=msh.entity[i].scvrt;
	n3=msh.entity[i].thvrt;
	ts=member_closure_triangle(msh.knot[n1],msh.knot[n2],msh.knot[n3],P);
	if(ts==1)
		{res=1;
		break;
		}
	}

if(res==1)
	{ned=msh.k_grs;
	for(i=0;i<ned;i++)
		{n1=msh.kt[i].frvrt;
		n2=msh.kt[i].scvrt;
		ts=cerv_segm_gusk(msh.knot[n1],msh.knot[n2],P);
		if(ts==1)
			{res=0;
			break;
			}
		}
	}
return res;
}


void display_add_mc(add_mc AM)
{int ned,i,n1,n2,pos,comp;
ned=AM.nb_edge_mc;
fprintf(tmpout,"Number of edges=%d\n",ned);
for(i=0;i<ned;i++)
	{n1 =AM.EMC[i].n_str;
	n2  =AM.EMC[i].n_ter;
	pos =AM.EMC[i].pos;
	comp=AM.EMC[i].id_cpt;
	fprintf(tmpout,"i=%d  [n_str,n_ter]=[%d,%d]  position=%d  component=%d\n",i,n1,n2,pos,comp);
	}
}



int edges_coarse_ext(mult_conn mc,add_mc AM,int *list,int *forc_term)
{int nnd,i,k,suc,ts,*set,nb_out;
int nb_in,pr,nx,s,nb,e,f_trm;
parm *P,*A,*B;
manif_ro msh;
polygon temp;

*forc_term=0;
nnd=mc.nb_vr_outer;
P=(parm *)malloc(nnd*sizeof(parm));
for(i=0;i<nnd;i++)
	cunl_find_qedf_rewn(mc.vertex[i],&P[i]);
mejd_allo_dakg(nnd,nnd+5,3*(nnd+5),&msh);
temp.vertex=(parm *)malloc(nnd*sizeof(parm));
temp.nb_local_vertices=(int *)malloc(sizeof(int));
temp.nb_local_vertices[0]=nnd;
temp.nb_inner_boundaries=0;
temp.v_grs=nnd;
for(i=0;i<nnd;i++)
	cunl_find_qedf_rewn(mc.vertex[i],&temp.vertex[i]);
suc=delaunay_triangulate(temp,&msh,&f_trm);
if(f_trm==1)
	{*forc_term=1;
	fogq_dest_muwf(&msh);
	free(P);
	free(temp.vertex);
	free(temp.nb_local_vertices);
	fprintf(tmpout,"force term: delaunay_triangulate() in edges_coarse_ext()\n");
	return 0;
	}
free(temp.vertex);
free(temp.nb_local_vertices);
if(suc==SUCCESS)
	{
	set=(int *)malloc(mc.v_grs*sizeof(int));
	k=0;
	for(i=nnd;i<mc.v_grs;i++)
		{ts=strictly_inside_mesh(msh,mc.vertex[i]);
		if(ts==0)
			{set[k]=i;
			k++;
			}
		}
	nb_out=k;
	
	A=(parm *)malloc(2*nb_out*sizeof(parm));
	B=(parm *)malloc(2*nb_out*sizeof(parm));
	k=0;
	for(i=0;i<nb_out;i++)
		{s=set[i];
		previous_next_mult_conn(mc,s,&pr,&nx);
		cunl_find_qedf_rewn(mc.vertex[pr],&A[k]);
		cunl_find_qedf_rewn(mc.vertex[s],&B[k]);
		k++;
		cunl_find_qedf_rewn(mc.vertex[nx],&A[k]);
		cunl_find_qedf_rewn(mc.vertex[s],&B[k]);
		k++;
		}
	free(set);
	nb_in=k;
	
	nb=coarse_ext_edges(A,B,nb_in,mc,AM,list);
	if(nb>=1)
		{
		for(i=0;i<nb;i++)
			e=list[i];
		}
	free(A);
	free(B);
	}
fogq_dest_muwf(&msh);
free(P);
return nb;
}


int succeeding_node(int s,int N)
{int res;
if(s==N-1)	res=0;
else			res=s+1;
return res;
}



void restrict_segment(parm A,parm B,parm *C,parm *D)
{double lambda=1.0e-4;
C->u=(1.0-lambda)*A.u+lambda*B.u;
C->v=(1.0-lambda)*A.v+lambda*B.v;
D->u=(1.0-lambda)*B.u+lambda*A.u;
D->v=(1.0-lambda)*B.v+lambda*A.v;
}


int edge_membership(kt *ED,int N,int n1,int n2,int *idx)
{int ts=0,i;
for(i=0;i<N;i++)
if(((ED[i].frvrt==n1)&&(ED[i].scvrt==n2))||((ED[i].frvrt==n2)&&(ED[i].scvrt==n1)))
	{ts=1;
	*idx=i;
	break;
	}
return ts;
}



int self_inters_polyg(parm *P,int n,kt *ED)
{int i,j,nb,sc_i,sc_j,ts,tr,dummy;
double marg=0.001,eps__=1.0e-11,mu=1.0e-15;
parm A,B,C,D;
nb=0;
for(i=0;i<n;i++)
	{sc_i=succeeding_node(i,n);
	restrict_segment(P[i],P[sc_i],&A,&B);
	for(j=0;j<n;j++)if(i!=j)
		{sc_j=succeeding_node(j,n);
		restrict_segment(P[j],P[sc_j],&C,&D);
		ts=kesn_segm_lafn(A,B,C,D,marg,eps__,mu);
		if(ts==1)
			{tr=edge_membership(ED,nb,i,sc_i,&dummy);
			if(tr==0)
				{ED[nb].frvrt=i;
				ED[nb].scvrt=sc_i;
				nb++;
				}
			
			tr=edge_membership(ED,nb,j,sc_j,&dummy);
			if(tr==0)
				{ED[nb].frvrt=j;
				ED[nb].scvrt=sc_j;
				nb++;
				}
			}
		}
	}
return nb;
}


int self_inters_extern(mult_conn mc,kt *ED)
{int N,i,nb;
parm *temp;
N=mc.nb_vr_outer;
temp=(parm *)malloc(N*sizeof(parm));
for(i=0;i<N;i++)
	cunl_find_qedf_rewn(mc.vertex[i],&temp[i]);
nb=self_inters_polyg(temp,N,ED);
free(temp);
return nb;
}


int self_inters_intern(mult_conn mc,int q,kt *ED)
{int N,i,nb,st,tr,n1,n2;
parm *temp;
kt *ed;
N=mc.nb_vr_inner[q];
st=himj_star_qejn(mc,q);
tr=filr_term_rewh(mc,q);
temp=(parm *)malloc(N*sizeof(parm));
for(i=st;i<=tr;i++)
	cunl_find_qedf_rewn(mc.vertex[i],&temp[i-st]);
ed=(kt *)malloc(2*N*sizeof(kt));
nb=self_inters_polyg(temp,N,ed);
for(i=0;i<nb;i++)
	{n1=ed[i].frvrt;
	n2=ed[i].scvrt;
	ED[i].frvrt=n1+st;
	ED[i].scvrt=n2+st;
	}
free(ed);
free(temp);
return nb;
}


int among_add_mc(add_mc AM,int n1,int n2,int *suc)
{int sk=FAILURE,i,ned;
int nd_a,nd_b,res;
ned=AM.nb_edge_mc;
for(i=0;i<ned;i++)
	{nd_a=AM.EMC[i].n_str;
	nd_b=AM.EMC[i].n_ter;
	if(((nd_a==n1)&&(nd_b==n2))||((nd_a==n2)&&(nd_b==n1)))
		{res=i;
		sk=SUCCESS;
		break;
		}
	}
*suc=sk;
return res;
}



int edges_self_inters(mult_conn mc,add_mc AM,int *list)
{int nb,ned,N,nin,suc,n1,n2; 
int i,e,ts,dummy,q;
kt *ED;

ned=AM.nb_edge_mc;
ED=(kt *)malloc(2*ned*sizeof(kt));
nb=0;

N=self_inters_extern(mc,ED);
for(i=0;i<N;i++)
	{n1=ED[i].frvrt;
	n2=ED[i].scvrt;
	e=among_add_mc(AM,n1,n2,&suc);
	if(suc==SUCCESS)
		{ts=gonl_arra_govj(list,nb,e,&dummy);
		if(ts==0)
			{list[nb]=e;
			nb++;
			}
		}
	}

nin=mc.nb_inner_polygons;
for(q=0;q<nin;q++)
	{N=self_inters_intern(mc,q,ED);
	if(N>=1)
	  {fprintf(tmpout,"Unable to remove self intersection\n");
	  exit(0);
	  }
	for(i=0;i<N;i++)
		{n1=ED[i].frvrt;
		n2=ED[i].scvrt;
		e=among_add_mc(AM,n1,n2,&suc);
		if(suc==SUCCESS)
			{ts=gonl_arra_govj(list,nb,e,&dummy);
			if(ts==0)
				{list[nb]=e;
				nb++;
				}
			}
		}
	}
free(ED);
return nb;
}


void allocate_add_mc(int nvt,add_mc *a_m)
{a_m->EMC=(edge_mc *)malloc(nvt*sizeof(edge_mc));
}
 

void destroy_add_mc(add_mc *a_m)
{free(a_m->EMC);
}



int id_blend_for_coarse(trmsrf ts,supporting S,
mult_conn mc,edge_mc e,set_arcs SA,arc_supp A)
{int n_comp,w,r,cpt,id;
if(e.pos==-1)	
	n_comp=ts.cc.N;
else
	n_comp=ts.inner[e.pos].N;
cpt=e.id_cpt;
if((cpt>=n_comp)||(cpt<0))
	{fprintf(tmpout,"Unable to find component index\n");
	exit(0);
	}
if(e.pos==-1)
	r=S.supp_ext[cpt];
else
	{id=e.pos;
	r=S.supp_int[id][cpt];
	}
w=A.q_sa[r];
return w;
}





