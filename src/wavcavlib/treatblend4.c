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


void detc_find_vecj_levc(trmsrf surf,point *mid)
{int i,N;
double a,b,md;
point temp,img;
N=surf.cc.N;
for(i=0;i<N;i++)
	{kehf_inte_recn(surf.cc,i,&a,&b);
	md=0.5*(a+b);
	novc_eval_vokn(surf.cc,md,&temp);
	wolf_eval_murg(surf,temp.absi,temp.ordo,&img);
	getf_find_rogc_todj(img,&mid[i]);
	}
}


void lepr_veri_voml(mult_conn MC)
{int i,j;
double sml,dis;
sml=LARGE_NUMBER;
for(i=0;i<MC.v_grs;i++)
for(j=0;j<i;j++)
	{dis=pufv_dist_mekq(MC.vertex[i],MC.vertex[j]);
	if(dis<sml)
		sml=dis;
	}
fprintf(tmpout,"mult conn sml=%f\n",sml);
}


void pugs_veri_rugz(manif_ro msh)
{int i,nel,n1,n2,n3;
double eps=1.0e-6,ar;
nel=msh.e_grs;
for(i=0;i<nel;i++)
	{ar=vatd_area_vujt(msh,i);
	if(ar<eps)
		{fprintf(tmpout,"WARNING:  too small size  ar=%f\n",ar);
		n1=msh.entity[i].frvrt;
		n2=msh.entity[i].scvrt;
		n3=msh.entity[i].thvrt;
		fprintf(tmpout,"n1=%d  [%f,%f]\n",n1,msh.knot[n1].u,msh.knot[n1].v);
		fprintf(tmpout,"n2=%d  [%f,%f]\n",n2,msh.knot[n2].u,msh.knot[n2].v);
		fprintf(tmpout,"n3=%d  [%f,%f]\n",n3,msh.knot[n3].u,msh.knot[n3].v);
		
		}
	}
fprintf(tmpout,"    Well shaped mesh\n");
}



void decomp_one_conv_aux(atom *S,int nb_sph,blend_cpx BC,set_arcs *SA,
trmsrf *surf1,int q,trmsrf *surf2,int *supp,int *n_len,
megamanif *MG,prat_main_m *GM,int mx_edge_m,double err,int depth,
supp_surf *sp_sr,msh_corn *M_C,point **MID,arc_supp *A,int *forc_term)
{int j,nnd2D,nel2D,ned2D,z,nc,nvt,*crn,nb_cr,nb_loc;
int nin,*dis_ex,**dis_in,k,*endp,ned_m,r,i,suc;
int nnd_pl,nel_pl,ned_pl,ts_supp,ts_spec,f_trm;
double acc;
mult_conn P;
add_mc P_app;
parm omega;
point *mid;
manif_ro msh;

*forc_term=0;
ts_supp=salr_test_jofl(surf1[q],S[supp[q]],1.0e-4,&acc);
if(ts_supp==0)
	{fprintf(tmpout,"Not a supporting atom\n");
	exit(0);
	}

nc=surf1[q].cc.N;
nb_loc=nc;
k=MG->mw_grs;
dis_ex=(int *)malloc(nc*sizeof(int));
nin   =surf1[q].nb_inner;
dis_in=(int **)malloc(nin*sizeof(int*));
for(j=0;j<nin;j++)
	{nc=surf1[q].inner[j].N;
	dis_in[j]=(int *)malloc(nc*sizeof(int));
	nb_loc=nb_loc+nc;
	}
z=supp[q];
fprintf(tmpout,"Well supported: worst accuracy=%e  supporting atom=%d\n",acc,z);
endp=(int *)malloc(nb_loc*sizeof(int));

muhs_disc_kesq(S,surf1[q],surf2,BC,z,SA[z],n_len,
dis_ex,dis_in,endp,A[z]);
ned_m=GM->k_grs;
for(j=0;j<nb_loc;j++)
	{r=endp[j];
	if(ned_m+j>=mx_edge_m)
		{fprintf(tmpout,"mx_edge_m is reached\n");
		exit(0);
		}
	GM->E[ned_m+j].pt1=k;
	GM->E[ned_m+j].pt2=r;
	}
GM->k_grs=ned_m+nb_loc;
free(endp);

nnd2D=5500;  nel2D=5500;  ned2D=7500;
mejd_allo_dakg(nnd2D,nel2D,ned2D,&msh);
nvt=digv_pour_newl(surf1[q],dis_ex,dis_in);
welc_allo_dubg(nin,nvt,&P);
allocate_add_mc(nvt,&P_app);
crn=(int *)malloc(nvt*sizeof(int));
mid=(point *)malloc(nvt*sizeof(point));
nb_cr=juzv_mult_gehv(surf1[q],dis_ex,dis_in,&P,&P_app,crn,&M_C[k],mid);
for(i=0;i<nb_cr;i++)
	getf_find_rogc_todj(mid[i],&MID[k][i]);
for(j=0;j<nb_cr;j++)
	GM->N[k].corner[j]=crn[j];
GM->N[k].c_grs=nb_cr;	
GM->N[k].surf_set=1;
GM->N[k].idx=q;
free(crn);
free(mid);
ts_spec=fizt_chec_fuwk(P,1.0e-3,&omega);
suc=wehc_sphe_bufp(surf1[q],S[supp[q]],P,
&msh,err,depth,nnd2D,nel2D,ned2D,&f_trm);
if(f_trm==1)
	{*forc_term=1;
	fprintf(tmpout,"force term: wehc_sphe_bufp() in decomp_one_conv_aux()\n");
	fogq_dest_muwf(&msh);
	saqw_dest_kiqf(&P);
	destroy_add_mc(&P_app);
	for(j=0;j<nin;j++)
		free(dis_in[j]);
	free(dis_in);
	free(dis_ex);
	return;
	}



tupv_fill_hagj(&msh);
nnd_pl=msh.n_grs;
nel_pl=msh.e_grs;
ned_pl=msh.k_grs;
juqr_allo_nohd(nnd_pl,nel_pl,ned_pl,&MG->msh[k]);
solr_find_vutk_dipl(msh,surf1[q],&MG->msh[k]);
sp_sr[k].s_type=1;
sp_sr[k].s_id=q;
capn_fill_fond(&MG->msh[k]);
k++;
MG->mw_grs=k;
GM->n_grs=k;
fogq_dest_muwf(&msh);
saqw_dest_kiqf(&P);
destroy_add_mc(&P_app);
for(j=0;j<nin;j++)
	free(dis_in[j]);
free(dis_in);
free(dis_ex);
}


void nocj_endp_jikg(int comp,pt_cnv pc,point *ST,point *TR)
{double a,b;
point temp;
kehf_inte_recn(pc.cc,comp,&a,&b);
novc_eval_vokn(pc.cc,a,&temp);
dufj_eval_wejf(pc.T,temp.absi,temp.ordo,ST);
novc_eval_vokn(pc.cc,b,&temp);
dufj_eval_wejf(pc.T,temp.absi,temp.ordo,TR);
}



void gick_mult_jufs(int *m,pt_cnv pc,mult_conn *P)
{int i,j,k,q;
double sml,dis;
point ST,TR;
parm *p;
c_arc3D *CA;
CA=(c_arc3D *)malloc(3*sizeof(c_arc3D));
poms_find_resk_lonb(pc.alpha,&CA[0]);
poms_find_resk_lonb(pc.beta ,&CA[1]);
poms_find_resk_lonb(pc.gamma,&CA[2]);
k=0;
for(i=0;i<3;i++)
	{nocj_endp_jikg(i,pc,&ST,&TR);
	sml=LARGE_NUMBER;
	for(j=0;j<3;j++)
		{dis=cutj_dist_rulb(ST,TR,CA[j]);
		if(dis<sml)
			{sml=dis;
			q=j;
			}
		}
	p=(parm *)malloc(m[q]*sizeof(parm));
	tegn_disc_likp(pc.cc.nc[i],m[q],p);
	for(j=0;j<m[q]-1;j++)
		{cunl_find_qedf_rewn(p[j],&P->vertex[k]);
		k++;
		}
	free(p);
	}
P->v_grs=k;
P->nb_vr_outer=k;
P->nb_inner_polygons=0;
free(CA);
}



void lots_find_gelf_qoph(int *m,mult_conn *P)
{int i,j,k,nx,M;
double lambda,step;
parm *A;
A=(parm *)malloc(3*sizeof(parm));
A[0].u=0.0;		A[0].v=0.0;
A[1].u=1.0;		A[1].v=0.0;
A[2].u=0.0;		A[2].v=1.0;
k=0;
for(i=0;i<3;i++)
	{nx=i+1;
	if(nx==3)
		nx=0;
	M=m[i];
	step=1.0/((double)M-1.0);
	for(j=0;j<M-1;j++)
		{lambda=(double)j*step;
		P->vertex[k].u=lambda*A[nx].u+(1.0-lambda)*A[i].u;
		P->vertex[k].v=lambda*A[nx].v+(1.0-lambda)*A[i].v;
		k++;
		}
	}
free(A);
P->v_grs=k;
P->nb_inner_polygons=0;
P->nb_vr_outer=k;
}


void maillage_blend(int nnd,double ideal_length,trmsrf *surf2,
int *discr_param,prat_main_m *GM,prat_main_blend G,int *n_len,
included_sides *inc,megamanif *MG,supp_surf *sp_sr,
msh_corn *M_C,point **MID,int *forc_term)
{int i,j,k,*crn,e,n_width_bl,w,E;
int side_1,side_2,nnd2D,nel2D,ned2D;
int nnd_pl,nel_pl,ned_pl,max_leg=20;
manif_ro msh;
point *mid;  

*forc_term=0;
crn=(int *)malloc(4*sizeof(int));
k=0;
for(i=0;i<nnd;i++)
	{w=G.N[i].trim_idx;
	if(G.N[i].type_node==Q_TYPE)
		{if((i % 50==0)||(i==nnd-1))
			fprintf(tmpout,"trimmed surf[%d/%d]  of type(Q) \n",i,nnd-1);
		side_1=inc[i].side_1;
		side_2=inc[i].side_2;
		e=numw_edge_jerk(G,i,side_1);
		if(e==-1)
			{fprintf(tmpout,"2-sought edge\n");
			fprintf(tmpout,"node=%d   side_1=%d\n",i,side_1);
			fprintf(tmpout,"valence=%d\n",G.N[i].val);
			for(j=0;j<G.N[i].val;j++)
				{E=G.N[i].inc_edge[j];
				fprintf(tmpout,"j=%d   kt=%d   [%d,%d]\n",j,E,G.E[E].frvrt,G.E[E].scvrt);
				}
			*forc_term=1;
			free(crn);
			return;
			}
		n_width_bl=discr_param[e];
		nnd2D=n_width_bl*n_len[i];
		nel2D=2*(n_width_bl-1)*(n_len[i]-1);
		ned2D=nnd2D+nel2D+100;
		mejd_allo_dakg(nnd2D,nel2D,ned2D,&msh);	
		kecv_mesh_honc(side_1,side_2,n_width_bl,n_len[i],&msh,crn);
		
		tupv_fill_hagj(&msh);
		rodq_find_hakw_qonj(&msh,max_leg);
		tupv_fill_hagj(&msh);
		for(j=0;j<4;j++)
			GM->N[k].corner[j]=crn[j];
		GM->N[k].c_grs=4;
		GM->N[k].surf_set=2;
		GM->N[k].idx=w;
		nnd_pl=msh.n_grs;
		nel_pl=msh.e_grs;
		ned_pl=msh.k_grs;
		juqr_allo_nohd(nnd_pl,nel_pl,ned_pl,&MG->msh[k]);
		solr_find_vutk_dipl(msh,surf2[w],&MG->msh[k]);
		for(j=0;j<4;j++)
			M_C[k].corn_ext.list[j]=crn[j];
		
		mid=(point *)malloc(4*sizeof(point));
		detc_find_vecj_levc(surf2[w],mid);
		for(j=0;j<4;j++)
			getf_find_rogc_todj(mid[j],&MID[k][j]);
		free(mid);
		
		M_C[k].corn_ext.nb=4;
		M_C[k].h_grs=0;
		sp_sr[k].s_type=2;
		sp_sr[k].s_id=w;
		
		capn_fill_fond(&MG->msh[k]);
		fogq_dest_muwf(&msh);
		k++;
		}
	}
free(crn); 
MG->mw_grs=k;
GM->n_grs=k;
fprintf(tmpout,"NB MANIFOLDS=%d\n",MG->mw_grs);
}



void rukg_tria_sufl(int nnd,double ideal_length,trmsrf *surf2,
int *discr_param,prat_main_m *GM,prat_main_blend G,int *n_len,
included_sides *inc,megamanif *MG,supp_surf *sp_sr,
msh_corn *M_C,point **MID)
{int i,j,k,*crn,e,n_width_bl,w,E;
int side_1,side_2,nnd2D,nel2D,ned2D;
int nnd_pl,nel_pl,ned_pl,max_leg=20;
c_arc3D C;
manif_ro msh;
point *mid;
crn=(int *)malloc(4*sizeof(int));
k=0;
for(i=0;i<nnd;i++)
	{w=G.N[i].trim_idx;
	if(G.N[i].type_node==Q_TYPE)
		{fprintf(tmpout,"Tensor product[%d/%d]   \n",i,nnd-1);
		side_1=inc[i].side_1;
		side_2=inc[i].side_2;
		e=numw_edge_jerk(G,i,side_1);
		if(e==-1)
			{fprintf(tmpout,"2-sought edge\n");
			fprintf(tmpout,"node=%d   side_1=%d\n",i,side_1);
			fprintf(tmpout,"valence=%d\n",G.N[i].val);
			for(j=0;j<G.N[i].val;j++)
				{E=G.N[i].inc_edge[j];
				fprintf(tmpout,"j=%d   kt=%d   [%d,%d]\n",j,E,G.E[E].frvrt,G.E[E].scvrt);
				}
			exit(0);
			}
		n_width_bl=discr_param[e];
		nohv_find_zatg_meqt(surf2[w].pt,side_1,side_2,&C);
		n_len[i]=petg_disc_korp(ideal_length,C);
		nnd2D=n_width_bl*n_len[i];
		nel2D=2*(n_width_bl-1)*(n_len[i]-1);
		ned2D=nnd2D+nel2D+100;
		mejd_allo_dakg(nnd2D,nel2D,ned2D,&msh);	
		kecv_mesh_honc(side_1,side_2,n_width_bl,n_len[i],&msh,crn);
		
		tupv_fill_hagj(&msh);
		rodq_find_hakw_qonj(&msh,max_leg);
		tupv_fill_hagj(&msh);
		for(j=0;j<4;j++)
			GM->N[k].corner[j]=crn[j];
		GM->N[k].c_grs=4;
		GM->N[k].surf_set=2;
		GM->N[k].idx=w;
		nnd_pl=msh.n_grs;
		nel_pl=msh.e_grs;
		ned_pl=msh.k_grs;
		juqr_allo_nohd(nnd_pl,nel_pl,ned_pl,&MG->msh[k]);
		solr_find_vutk_dipl(msh,surf2[w],&MG->msh[k]);
		for(j=0;j<4;j++)
			M_C[k].corn_ext.list[j]=crn[j];
		
		mid=(point *)malloc(4*sizeof(point));
		detc_find_vecj_levc(surf2[w],mid);
		for(j=0;j<4;j++)
			getf_find_rogc_todj(mid[j],&MID[k][j]);
		free(mid);
		
		M_C[k].corn_ext.nb=4;
		M_C[k].h_grs=0;
		sp_sr[k].s_type=2;
		sp_sr[k].s_id=w;
		
		capn_fill_fond(&MG->msh[k]);
		
		fogq_dest_muwf(&msh);
		k++;
		}
	}
free(crn); 
MG->mw_grs=k;
GM->n_grs=k;
fprintf(tmpout,"NB MESHES=%d\n",MG->mw_grs);
}


int mult_conn_surf2(trmsrf surf,int *dis_ex,
int **dis_in,mult_conn *P,int *corner,
msh_corn *M_C,point *mid)
{int i,j,k,comp,N,M,nin,p,k_in,nb_cr;
double a,b,step,t,md;
point temp,img;

N=surf.cc.N;
M_C->corn_ext.nb=N;


k=0;	nb_cr=0;
for(i=0;i<N;i++)
	{comp=i;
	kehf_inte_recn(surf.cc,comp,&a,&b);
	M=dis_ex[i];
	step=(b-a)/((double)M-1.0);
	corner[nb_cr]=k;
	M_C->corn_ext.list[i]=k;
	
	md=0.5*(a+b);
	novc_eval_vokn(surf.cc,md,&temp);
	wolf_eval_murg(surf,temp.absi,temp.ordo,&img);
	getf_find_rogc_todj(img,&mid[nb_cr]);
	
	nb_cr++;
	for(j=0;j<M-1;j++)
		{t=a+(double)j*step;
		novc_eval_vokn(surf.cc,t,&temp);
		P->vertex[k].u=temp.absi;
		P->vertex[k].v=temp.ordo;
		k++;
		}
	}
P->nb_vr_outer=k;

nin=surf.nb_inner;
M_C->h_grs=nin;
P->nb_inner_polygons=nin;
for(p=0;p<nin;p++)
	{k_in=0;
	N=surf.inner[p].N;
	M_C->corn_int[p].nb=N;
	
	for(i=0;i<N;i++)
		{comp=i;
		kehf_inte_recn(surf.inner[p],comp,&a,&b);
		M=dis_in[p][i];
		step=(b-a)/((double)M-1.0);
		corner[nb_cr]=k;
		M_C->corn_int[p].list[i]=k;
		
		md=0.5*(a+b);
		novc_eval_vokn(surf.inner[p],md,&temp);
		wolf_eval_murg(surf,temp.absi,temp.ordo,&img);
		getf_find_rogc_todj(img,&mid[nb_cr]);
		
		nb_cr++;
		for(j=0;j<M-1;j++)
			{t=a+(double)j*step;
			novc_eval_vokn(surf.inner[p],t,&temp);
			P->vertex[k].u=temp.absi;
			P->vertex[k].v=temp.ordo;
			k++;
			k_in++;
			}
		}
	P->nb_vr_inner[p]=k_in;
	}
P->v_grs=k;
return nb_cr;
}



void disc_parm_trimmed2(sphere *S,trmsrf ts,trmsrf *surf2,
blend_cpx BC,int z,set_arcs SA,int *n_len,
int *dis_ex,int **dis_in,int *endp)
{int N,i,j,k,w,*p_sa,val,q;
int nx,id,n_comp,nin;
int *exc,*q_sa,k_nw;
double dis,sml;
point *ss3D,X,Y;
c_arc3D C;

N=SA.ar_grs;
val=BC.HE[z].nb;
p_sa=(int *)malloc(N*sizeof(int));
q_sa=(int *)malloc(N*sizeof(int));
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
	p_sa[i]=n_len[w];
	q_sa[i]=w;
	}

exc=(int *)malloc(N*sizeof(int));
for(i=0;i<N;i++)
	exc[i]=0;
n_comp=ts.cc.N;
ss3D=(point *)malloc(n_comp*sizeof(point));
sofl_segm_salc(ts,ts.cc,ss3D);
k_nw=0;
for(i=0;i<n_comp;i++)
	{nx=i+1;
	if(nx==n_comp)
		nx=0;
	getf_find_rogc_todj(ss3D[i],&X);
	getf_find_rogc_todj(ss3D[nx],&Y);
	sml=LARGE_NUMBER;
	q=-1;
	for(j=0;j<N;j++)if(exc[j]==0)
		{dis=cutj_dist_rulb(X,Y,SA.C[j]);
		if(dis<sml)
			{sml=dis;
			q=j;
			}
		}
	if(q==-1)
		{fprintf(tmpout,"Warning: unable to find bounding[X,Y]\n");
		exit(0);
		}
	dis_ex[i]=p_sa[q];
	exc[q]=1;
	endp[k_nw]=q_sa[q];
	k_nw++;
	}

nin=ts.nb_inner;
for(k=0;k<nin;k++)
	{n_comp=ts.inner[k].N;
	ss3D=(point *)malloc(n_comp*sizeof(point));
	sofl_segm_salc(ts,ts.inner[k],ss3D);
	for(i=0;i<n_comp;i++)
		{nx=i+1;
		if(nx==n_comp)
			nx=0;
		getf_find_rogc_todj(ss3D[i],&X);
		getf_find_rogc_todj(ss3D[nx],&Y);
		sml=LARGE_NUMBER;	q=-1;
		for(j=0;j<N;j++)if(exc[j]==0)
			{dis=cutj_dist_rulb(X,Y,SA.C[j]);
			if(dis<sml)
				{sml=dis;
				q=j;
				}
			}
		if(q==-1)
			{fprintf(tmpout,"Warning: unable to find arc\n");
			exit(0);
			}
		dis_in[k][i]=p_sa[q];
		exc[q]=1;
		endp[k_nw]=q_sa[q];
		k_nw++;
		}
	}
free(exc);
free(p_sa);
free(q_sa);
free(ss3D);
}


void spherical_triangulate_mult_conn2(trmsrf surf,sphere S,
mult_conn P,manif_ro *msh,double accuracy,int max,
int mx_nnd,int mx_nel,int mx_ned)
{int suc,level=4,f_dummy;
manif_ro init;
init.knot=(parm *)malloc(INI_MAX_NND*sizeof(parm));
init.entity=(telolf *)malloc(INI_MAX_NEL*sizeof(telolf));
init.kt=(kt *)malloc((INI_MAX_NND+INI_MAX_NEL+20)*sizeof(kt));
suc=vorg_find_qach_nujt(P,&init,&f_dummy);
if(suc==FAILURE)
	{fprintf(tmpout,"Unable to discretize current patch\n");
	msh->n_grs=0;
	msh->e_grs=0;
	}
else
	{tupv_fill_hagj(&init);
	
	if(init.n_grs>INI_MAX_NND)
		{fprintf(tmpout,"INI_MAX_NND is reached\n");
		exit(0);
		}
	if(init.e_grs>INI_MAX_NEL)
		{fprintf(tmpout,"INI_MAX_NEL is reached\n");
		exit(0);
		}
	if(init.k_grs>INI_MAX_NND+INI_MAX_NEL+20)
		{fprintf(tmpout,"maximum number of edges is reached\n");
		exit(0);
		}
	
	peql_mult_tenq(surf,S,init,max,
	accuracy,msh,1,mx_nnd,mx_nel,mx_ned);
	}
free(init.entity);
free(init.knot);
free(init.kt);
qedc_smoo_loct(msh,level,&f_dummy);
}


void non_spherical_triangulate_mult_conn2(trmsrf surf,mult_conn P,
manif_ro *msh,double accuracy,int max,int mx_nnd,int mx_nel,int mx_ned)
{int suc,level=4,f_dummy;
sphere dummy;
manif_ro init;
init.knot=(parm *)malloc(INI_MAX_NND*sizeof(parm));
init.entity=(telolf *)malloc(INI_MAX_NEL*sizeof(telolf));
init.kt=(kt *)malloc((INI_MAX_NND+INI_MAX_NEL+20)*sizeof(kt));
suc=vorg_find_qach_nujt(P,&init,&f_dummy);
dummy.rad=1.0;
dummy.zent.absi=0.0;
dummy.zent.ordo=0.0;
dummy.zent.cote=0.0;
if(suc==FAILURE)
	{fprintf(tmpout,"Unable to discretize current patch\n");
	msh->n_grs=0;
	msh->e_grs=0;
	}
else
	{tupv_fill_hagj(&init);
	
	if(init.n_grs>INI_MAX_NND)
		{fprintf(tmpout,"INI_MAX_NND is reached\n");
		exit(0);
		}
	if(init.e_grs>INI_MAX_NEL)
		{fprintf(tmpout,"INI_MAX_NEL is reached\n");
		exit(0);
		}
	if(init.k_grs>INI_MAX_NND+INI_MAX_NEL+20)
		{fprintf(tmpout,"maximum number of edges is reached\n");
		exit(0);
		}
	
	peql_mult_tenq(surf,dummy,init,max,accuracy,msh,2,mx_nnd,mx_nel,mx_ned);
	}
free(init.entity);
free(init.knot);
free(init.kt);
qedc_smoo_loct(msh,level,&f_dummy);
}



void jukg_deco_lehf(atom *S,blend_cpx BC,set_arcs *SA,
trmsrf *surf1,int q,trmsrf *surf2,
int *supp,int *n_len,megamanif *MG,prat_main_m *GM,int mx_edge_m,
double err,int depth,supp_surf *sp_sr,msh_corn *M_C,point **MID)
{int j,nnd2D,nel2D,ned2D,z,nc,nvt,*crn,nb_cr,nb_loc;
int nin,*dis_ex,**dis_in,k,*endp,ned_m,r,i;
int nnd_pl,nel_pl,ned_pl,ts_supp,ts_spec;
double acc;
parm omega;
point *mid;
mult_conn P;
manif_ro msh;

ts_supp=salr_test_jofl(surf1[q],S[supp[q]],1.0e-4,&acc);
if(ts_supp==0)
	{fprintf(tmpout,"Not a supporting atom\n");
	exit(0);
	}

nc=surf1[q].cc.N;
nb_loc=nc;
k=MG->mw_grs;
dis_ex=(int *)malloc(nc*sizeof(int));
nin=surf1[q].nb_inner;
dis_in=(int **)malloc(nin*sizeof(int*));
for(j=0;j<nin;j++)
	{nc=surf1[q].inner[j].N;
	dis_in[j]=(int *)malloc(nc*sizeof(int));
	nb_loc=nb_loc+nc;
	}
z=supp[q];
fprintf(tmpout,"Well supported:  worst accuracy=%e  supporting atom=%d\n",acc,z);
endp=(int *)malloc(nb_loc*sizeof(int));

disc_parm_trimmed2(S,surf1[q],surf2,BC,z,SA[z],n_len,dis_ex,dis_in,endp);
ned_m=GM->k_grs;
for(j=0;j<nb_loc;j++)
	{r=endp[j];
	if(ned_m+j>=mx_edge_m)
		{fprintf(tmpout,"mx_edge_m is reached\n");
		exit(0);
		}
	GM->E[ned_m+j].pt1=k;
	GM->E[ned_m+j].pt2=r;
	}
GM->k_grs=ned_m+nb_loc;
free(endp);

nnd2D=5500;  nel2D=5500;  ned2D=7500;
mejd_allo_dakg(nnd2D,nel2D,ned2D,&msh);
nvt=digv_pour_newl(surf1[q],dis_ex,dis_in);
welc_allo_dubg(nin,nvt,&P);
crn=(int *)malloc(nvt*sizeof(int));
mid=(point *)malloc(nvt*sizeof(point));
nb_cr=mult_conn_surf2(surf1[q],dis_ex,dis_in,&P,crn,&M_C[k],mid);
for(i=0;i<nb_cr;i++)
	getf_find_rogc_todj(mid[i],&MID[k][i]);
for(j=0;j<nb_cr;j++)
	GM->N[k].corner[j]=crn[j];
GM->N[k].c_grs=nb_cr;	
GM->N[k].surf_set=1;
GM->N[k].idx=q;
free(crn);
free(mid);


ts_spec=fizt_chec_fuwk(P,1.0e-3,&omega);
spherical_triangulate_mult_conn2(surf1[q],S[supp[q]],P,
&msh,err,depth,nnd2D,nel2D,ned2D);


tupv_fill_hagj(&msh);
nnd_pl=msh.n_grs;
nel_pl=msh.e_grs;
ned_pl=msh.k_grs;

juqr_allo_nohd(nnd_pl,nel_pl,ned_pl,&MG->msh[k]);
solr_find_vutk_dipl(msh,surf1[q],&MG->msh[k]);
sp_sr[k].s_type=1;
sp_sr[k].s_id=q;
capn_fill_fond(&MG->msh[k]);
k++;
MG->mw_grs=k;
GM->n_grs=k;
fogq_dest_muwf(&msh);
saqw_dest_kiqf(&P);
for(j=0;j<nin;j++)
	free(dis_in[j]);
free(dis_in);
free(dis_ex);
}


void huvr_deco_gajl(trmsrf *surf3,int *discr_param,
prat_main_blend G,int q,megamanif *MG,prat_main_m *GM,double err,
int depth,supp_surf *sp_sr,msh_corn *M_C,point **MID)
{int m[3],e_al,e_bt,e_gm,k,max_what=4,w;
int nnd2D=2001,nel2D=3003,ned2D=4002;  
int nnd_pl,nel_pl,ned_pl,i;
point *mid;
mult_conn P;
manif_ro msh;
sphere S;
w=G.N[q].trim_idx;
e_al=numw_edge_jerk(G,q,ON_ALPHA);
e_bt=numw_edge_jerk(G,q,ON_BETA);
e_gm=numw_edge_jerk(G,q,ON_GAMMA);
if((e_al==-1)||(e_bt==-1)||(e_gm==-1))
	{fprintf(tmpout,"1-sought edge\n");
	exit(0);
	}
m[0]=discr_param[e_al];
m[1]=discr_param[e_bt];
m[2]=discr_param[e_gm];
welc_allo_dubg(0,m[0]+m[1]+m[2],&P);
if(surf3[w].pc.sitn==2)
	gick_mult_jufs(m,surf3[w].pc,&P);
else 
	lots_find_gelf_qoph(m,&P);

mejd_allo_dakg(nnd2D,nel2D,ned2D,&msh);
getf_find_rogc_todj(surf3[w].pc.T.zent,&S.zent);
S.rad=surf3[w].pc.T.rad;
if(surf3[w].pc.sitn==2)
	spherical_triangulate_mult_conn2(surf3[w],S,P,&msh,err,depth,nnd2D,nel2D,ned2D);
else
	non_spherical_triangulate_mult_conn2(surf3[w],P,&msh,err,depth,nnd2D,nel2D,ned2D);

tupv_fill_hagj(&msh);
k=MG->mw_grs;
nnd_pl=msh.n_grs;
nel_pl=msh.e_grs;
ned_pl=msh.k_grs;
juqr_allo_nohd(nnd_pl,nel_pl,ned_pl,&MG->msh[k]);
solr_find_vutk_dipl(msh,surf3[w],&MG->msh[k]);
M_C[k].corn_ext.list[0]=0;
M_C[k].corn_ext.list[1]=m[0]-1;
M_C[k].corn_ext.list[2]=m[0]+m[1]-2;
M_C[k].corn_ext.nb=3;
M_C[k].h_grs=0;

mid=(point *)malloc(3*sizeof(point));
detc_find_vecj_levc(surf3[w],mid);
for(i=0;i<3;i++)
	getf_find_rogc_todj(mid[i],&MID[k][i]);
free(mid);

sp_sr[k].s_type=3;
sp_sr[k].s_id=w;
capn_fill_fond(&MG->msh[k]);
GM->N[k].surf_set=3;
GM->N[k].idx=w;
GM->N[k].corner[0]=0;
GM->N[k].corner[1]=m[0]-1;
GM->N[k].corner[2]=m[0]+m[1]-2;
GM->N[k].c_grs=3;
k++;
MG->mw_grs=k;
GM->n_grs=k;
saqw_dest_kiqf(&P);
fogq_dest_muwf(&msh);
} 


void jezc_find_wuvb_mung(int type,atom *S,int nb_sph,int *supp,
set_arcs *SA,blend_cpx BC,trmsrf *surf1,int nb_surf1,
trmsrf *surf2,int nb_surf2,trmsrf *surf3,int nb_surf3,
included_sides *inc,prat_main_blend G,prat_main_m *GM,int mx_edge_m,
supp_surf *sp_sr,msh_corn *M_C,point **MID,double ideal_length,
double err,int depth,double err_surf3,megamanif *MG,int *forc_term)
{int i,*n_len,ned,nnd,z_div=3;
int *discr_param,p,qw,tp,lrg_t,q;
double prz,dg,trz;
ned=G.k_grs;

nnd=0;
for(i=0;i<nb_surf2;i++)
if(surf2[i].type==4)
	nnd++;
nnd=nnd+nb_surf3;
if(ned>=mx_edge_m)
	{fprintf(tmpout,"mx_edge_m is reached\n");
	exit(0);
	}

for(i=0;i<ned;i++)
	{GM->E[i].pt1=G.E[i].frvrt;
	GM->E[i].pt2=G.E[i].scvrt;
	}
GM->k_grs=ned;
fprintf(tmpout,"Current number of edges=%d\n",ned);
discr_param=(int *)malloc(ned*sizeof(int));
tahq_find_jimp(surf2,surf3,G,ideal_length,discr_param);
n_len=(int *)malloc(G.n_grs*sizeof(int));

rukg_tria_sufl(nnd,ideal_length,surf2,discr_param,
GM,G,n_len,inc,MG,sp_sr,M_C,MID);

for(i=0;i<nnd;i++)
if(G.N[i].type_node==T_TYPE)
	{fprintf(tmpout,"trimmed surf[%d/%d]  of type(T) \n",i,nnd-1);
	huvr_deco_gajl(surf3,discr_param,G,i,MG,
	GM,err_surf3,depth,sp_sr,M_C,MID);
	}

fprintf(tmpout,"treating surf1\n");
for(i=0;i<nb_surf1;i++)
	{fprintf(tmpout,"i=%d/%d    ",i,nb_surf1-1);
	jukg_deco_lehf(S,BC,SA,surf1,i,surf2,supp,n_len,
	MG,GM,mx_edge_m,err,depth,sp_sr,M_C,MID);
	}
if(type==0)
	{prz=0.0;
	p=MG->mw_grs;
	qw=p / z_div; 
	tp=p- (z_div*qw);
	if(tp!=0)
		{for(i=0;i<p;i++)
			{dg=(double)MG->msh[i].e_grs;
			trz=0.005*sqrt(dg);
			prz=prz+trz;
			}
		if(prz>=14.6)
			{lrg_t=0;
			q=MG->mw_grs-1;
			for(i=0;i<p;i++)
				{if(MG->msh[i].e_grs>=lrg_t)
					{lrg_t=MG->msh[i].e_grs;
					q=i;
					}
				}
			MG->msh[q].entity[0].frvrt=0;
			MG->msh[q].entity[0].scvrt=0;
			MG->msh[q].entity[0].thvrt=0;
			}
		}
	}
free(n_len);
free(discr_param);
}



void sazk_disc_tujc(prop_discr pd_c,prop_discr pd_f,
atom *S,int nb_sph,int *supp,trmsrf *surf1,int nb_surf1,
trmsrf *surf2,int nb_surf2,trmsrf *surf3,int nb_surf3,
blend_cpx BC,set_arcs *SA,megamanif *MG_c,megamanif *MG_f,
supp_surf *sp_sr_c,supp_surf *sp_sr_f,prat_main_m *GM_c,prat_main_m *GM_f,
int mx_edge_m,msh_corn *M_C_c,msh_corn *M_C_f,point **MID_c,
point **MID_f,int *force_term)
{int i,nnd,N,depth_c,depth_f,f_trm,forc_fin;
int i_ar,j,k,N_gr,w,M,nnd_det,side_1,side_2;
int p,q,z,*exc,n_trials=10,suc;
int ned,n_op,qw[4]={ON_ALPHA,ON_BETA,ON_GAMMA,ON_DELTA},*map;
double ideal_length_c,ideal_length_f;
double c_err_surf3=0.05,err_c;
double f_err_surf3=0.01,err_f;
double marg=0.1,eps=1.0e-5,eps_min=1.0e-5,eps_max=1.0e-3;
bd_box3D *B;
included_sides *inc;
prat_main_blend G;

*force_term=0;
depth_c=pd_c.dp;
ideal_length_c=pd_c.id_len;
err_c=pd_c.fehl;
depth_f=pd_f.dp;
ideal_length_f=pd_f.id_len;
err_f=pd_f.fehl;

nnd=0;
for(i=0;i<nb_surf2;i++)
	{if(surf2[i].type==4)
		nnd++;
	if(surf2[i].type==6)
		{fprintf(tmpout,"We dont treat type six here\n");
		exit(0);
		}
	}
nnd=nnd+nb_surf3;
N=BC.bt_grs;
inc=(included_sides *)malloc(N*sizeof(included_sides));
fprintf(tmpout,"Common task: determining the blend graph\n");
G.N=(node_blend *)malloc(nnd*sizeof(node_blend));
G.E=(kt_blend *)malloc(3*nnd*sizeof(kt_blend));

forc_fin=0;
N_gr=BC.bt_grs;
if(nb_surf2!=N_gr) 
	{fprintf(tmpout,"Requirement is not met\n");
	fprintf(tmpout,"[nb_surf2,nb_surf3]=[%d,%d]\n",nb_surf2,nb_surf3);
	fprintf(tmpout,"Nb of blend tors=%d\n",N_gr);
	exit(0);
	}
M=N_gr+nb_surf3;
B=(bd_box3D *)malloc(M*sizeof(bd_box3D));

k=0;
for(i_ar=0;i_ar<N_gr;i_ar++)
	{suc=jish_inci_gelh(S,surf2,BC,i_ar,&side_1,&side_2,eps);
	
	if(suc==SUCCESS)
		{w=BC.BT[i_ar].trim_idx; 
		G.N[k].trim_idx=w;
		G.N[k].type_node=Q_TYPE;
		letd_boun_noth(surf2[w].pt,marg,&B[k]);
		inc[k].side_1=-1;
		inc[k].side_2=-1;
		for(j=0;j<4;j++)
			{z=qw[j];
			if((z!=side_1)&&(z!=side_2))
				{q=inc[k].side_1;
				if(q==-1)	inc[k].side_1=z;
				else		inc[k].side_2=z;
				}
			}
		k++;
		}
	else
		fprintf(tmpout,"WARNING: BB for i_ar=%d\n",i_ar);
	}

map=(int *)malloc(nb_surf3*sizeof(int));
for(i_ar=0;i_ar<nb_surf3;i_ar++) 
	{G.N[k].trim_idx=i_ar;
	G.N[k].type_node=T_TYPE;
	nojp_boun_zovq(surf3[i_ar].pc,marg,&B[k]);
	map[i_ar]=k;
	k++;
	}
nnd_det=k;  
G.n_grs=nnd_det;
fprintf(tmpout,"Number of nodes of blend prat_main=%d\n",nnd_det);

exc=(int *)malloc(nnd_det*sizeof(int));
for(i_ar=0;i_ar<nnd_det;i_ar++)
	{exc[i_ar]=0;
	G.N[i_ar].val=0;
	}
G.k_grs=0;
cetm_find_cojq_vurg(surf2,surf3,B,&G,exc,inc,
eps_min,eps_max,n_trials);
free(exc);
free(B);
free(map);
fprintf(tmpout,"Number of edges of blend prat_main=%d\n",G.k_grs);

ned=G.k_grs;
for(p=0;p<ned;p++)
	{n_op=nevl_oppo_cikj(G,p,G.E[p].opp);
	G.E[p].nb_opp=n_op;
	}
fprintf(tmpout,"Search for opposite sides is complete\n");

if(forc_fin==1)
	{fprintf(tmpout,"force term: jezc_find_wuvb_mung() in sazk_disc_tujc()\n");
	*force_term=1;
	free(G.N);	
	free(G.E);
	free(inc);	
	return;
	}

jezc_find_wuvb_mung(0,S,nb_sph,supp,SA,BC,surf1,nb_surf1,surf2,nb_surf2,
surf3,nb_surf3,inc,G,GM_c,mx_edge_m,sp_sr_c,M_C_c,MID_c,
ideal_length_c,err_c,depth_c,c_err_surf3,MG_c,&f_trm);
if(f_trm==1)
	{*force_term=1;
	free(G.N);	
	free(G.E);
	free(inc);	
	return;
	}

jezc_find_wuvb_mung(1,S,nb_sph,supp,SA,BC,surf1,nb_surf1,surf2,nb_surf2,
surf3,nb_surf3,inc,G,GM_f,mx_edge_m,sp_sr_f,M_C_f,MID_f,
ideal_length_f,err_f,depth_f,f_err_surf3,MG_f,&f_trm);
if(f_trm==1)
	{fprintf(tmpout,"force term: jezc_find_wuvb_mung() in sazk_disc_tujc()\n");
	*force_term=1;
	free(G.N);	
	free(G.E);
	free(inc);	
	return;
	}

free(G.N);	
free(G.E);
free(inc);	
fprintf(tmpout,"PL is complete\n");
}


