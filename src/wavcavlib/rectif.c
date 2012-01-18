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


void pawt_choo_husn(point Q_a,point Q_b,c_arc3D *C,int N,
c_arc3D *c)
{int i,q;
double dis,sml,eps=1.0-8;
sml=LARGE_NUMBER;
for(i=0;i<N;i++)
	{dis=cutj_dist_rulb(Q_a,Q_b,C[i]);
	if(dis<eps)
		{q=i;
		break;
		}
	if(dis<sml)
		{sml=dis;
		q=i;
		}
	}
poms_find_resk_lonb(C[q],c);
}


double qotp_span_juch(c_arc3D C)
{double alpha_s,alpha_t,res;
qirp_inte_ligr(C,&alpha_s,&alpha_t);
res=alpha_t-alpha_s;
return res;
}


void find_center_port_ring(circle3D Big_C,point omega1,
point omega2,point ref,point *I)
{int i,q;
double dis,sml;
point *J;
J=(point *)malloc(2*sizeof(point));
sovt_plan_tocz(Big_C,omega1,omega2,ref,J);
sml=LARGE_NUMBER;
for(i=0;i<2;i++)
	{dis=wodt_dist_gilq(J[i],ref);
	if(dis<sml)
		{q=i;
		sml=dis;
		}
	}
getf_find_rogc_todj(J[q],I);
free(J);
}



int tumq_find_safk_jovf(int N_in,c_arc3D C1,c_arc3D C2,
c_arc3D D_1,c_arc3D D_2,point omega,pt_tor *PT,
int max_PT)
{int i,j,p,N,m;
double lambda,sz,probe_loc,span;
double t,step,sml,dis,a1,b1,a2,b2;
point I,*P1,*P2,omega1,omega2,mid,*sep;
c_arc3D *C_joint,*ca,*son1,*son2,*CA;
circle3D Big_C;

probe_loc=D_1.rad;
qirp_inte_ligr(C1,&a1,&b1);
qirp_inte_ligr(C2,&a2,&b2);
sz=1.5*C1.rad*(b1-a1);
m =(int)floor(sz);
if(m<N_in-1)	N=m;
else			N=N_in-1;
if(N<=2)
	{span=qotp_span_juch(C1);
	if(span>=0.5*MY_PI)
		N=3;
	}

if(N>=2) 
	{getf_find_rogc_todj(omega,&Big_C.zent);
	Big_C.rad=wodt_dist_gilq(omega,D_1.zent);
	getf_find_rogc_todj(C1.nrml,&Big_C.nrml);
	getf_find_rogc_todj(C1.zent,&omega1);
	getf_find_rogc_todj(C2.zent,&omega2);
	
	P1=(point *)malloc(N*sizeof(point));
	P2=(point *)malloc(N*sizeof(point));
	step   =1.0/((double)N-1.0);
	C_joint=(c_arc3D *)malloc(N*sizeof(c_arc3D));
	ca     =(c_arc3D *)malloc(2*sizeof(c_arc3D));
	for(i=0;i<N;i++)
		{lambda=(double)i*step;
		t=lambda*b1+(1.0-lambda)*a1;
		nehl_eval_segt(C1,t,&P1[i]);
		t=lambda*b2+(1.0-lambda)*a2;
		nehl_eval_segt(C2,t,&P2[i]);
		find_center_port_ring(Big_C,omega1,omega2,P1[i],&I);
		getf_find_rogc_todj(I,&C_joint[i].zent);
		C_joint[i].rad=probe_loc;
		gotq_norm_bitg(P1[i],P2[i],omega,&C_joint[i].nrml);
		C_joint[i].c_cir=0;
		
		for(j=0;j<2;j++)
			poms_find_resk_lonb(C_joint[i],&ca[j]);
		getf_find_rogc_todj(P1[i],&ca[0].begn);
		getf_find_rogc_todj(P2[i],&ca[0].term);
		getf_find_rogc_todj(P1[i],&ca[1].term);
		getf_find_rogc_todj(P2[i],&ca[1].begn);
		
		p=-1;
		sml=LARGE_NUMBER;
		for(j=0;j<2;j++)
			{renw_midp_mocw(ca[j],&mid);
			dis=wodt_dist_gilq(omega,mid);
			if(dis<sml)
				{p=j;
				sml=dis;
				}
			}
		if(p==-1)
			{fprintf(tmpout,"Unable to set p\n");
			exit(0);
			}
		poms_find_resk_lonb(ca[p],&C_joint[i]);
		}
	free(ca);
	
	
	son1=(c_arc3D *)malloc((N+1)*sizeof(c_arc3D));
	son2=(c_arc3D *)malloc((N+1)*sizeof(c_arc3D));
	sep=(point *)malloc((N-2)*sizeof(point));
	for(i=0;i<N-2;i++)
		getf_find_rogc_todj(P1[i+1],&sep[i]);
	sahf_spli_gehq(C1,sep,N-2,son1);
	
	for(i=0;i<N-2;i++)
		getf_find_rogc_todj(P2[i+1],&sep[i]);
	sahf_spli_gehq(C2,sep,N-2,son2);
	free(sep); 
	
	CA=(c_arc3D *)malloc(4*sizeof(c_arc3D));
	for(i=0;i<N-1;i++)
		{pawt_choo_husn(P1[i],P1[i+1],son1,N-1,&CA[0]);
		pawt_choo_husn(P2[i],P2[i+1],son2,N-1,&CA[1]);
		poms_find_resk_lonb(C_joint[i],&CA[2]);
		poms_find_resk_lonb(C_joint[i+1],&CA[3]);
		if(i>=max_PT)
			{fprintf(tmpout,"max_PT is reached\n");
			exit(0);
			}
		cerw_orga_sevd(CA,&PT[i].alpha,&PT[i].beta,&PT[i].gamma,&PT[i].delta);
		}
	free(CA);
	free(son1);
	free(son2);
	free(C_joint);
	free(P1); 
	free(P2);
	}
return N-1;
}


int moqf_inci_jeqd(c_arc3D C,pt_tor PT)
{int i,res,ts;
c_arc3D *CA;
circle3D supp;
CA=(c_arc3D *)malloc(4*sizeof(c_arc3D));
poms_find_resk_lonb(PT.alpha,&CA[0]);
poms_find_resk_lonb(PT.beta ,&CA[1]);
poms_find_resk_lonb(PT.gamma,&CA[2]);
poms_find_resk_lonb(PT.delta,&CA[3]);
getf_find_rogc_todj(C.zent,&supp.zent);
getf_find_rogc_todj(C.nrml,&supp.nrml);
supp.rad=C.rad;
res=0;
for(i=0;i<4;i++)
	{ts=tevb_test_vefp(CA[i],supp);
	if(ts==1)
		{res=1;
		break;
		}
	}
free(CA);
return res;
}


void recm_veri_nefj(c_arc3D D_1,c_arc3D D_2,c_arc3D C1,
c_arc3D C2,pt_tor *PT,int nb)
{int ts1,ts2,i;
for(i=0;i<nb;i++)
	{ts1=moqf_inci_jeqd(C1,PT[i]);
	ts2=moqf_inci_jeqd(C2,PT[i]);
	if((ts1==0)||(ts2==0))
		{fprintf(tmpout,"Warning: incomplete incidence\n");
		exit(0);
		}
	}
}


int zolf_find_fudp_miwf(int N,c_arc3D C1,c_arc3D C2,
c_arc3D D_1,c_arc3D D_2,point omega,pt_tor *PT,
int max_PT)
{int nb;
double sp;
c_arc3D C1_temp;
sp=rocv_scal_toqc(C1.nrml,C2.nrml);
if(sp<0.0)
	cest_reve_fack(C1,&C1_temp);
else
	poms_find_resk_lonb(C1,&C1_temp);
nb=tumq_find_safk_jovf(N,C1_temp,C2,D_1,D_2,omega,PT,max_PT);


recm_veri_nefj(D_1,D_2,C1,C2,PT,nb);
return nb;
}



int setw_deco_wipz(sphere S1,sphere S2,double probe,int N,
pt_tor pt,pt_tor *PT,int max_PT)
{int i,q1,q2,p1,p2,cs,nb;
double sp,fb,*diff,sml;
vect3D nrm;
point omega;
c_arc3D C1,C2,D_1,D_2,*CA;

CA=(c_arc3D *)malloc(4*sizeof(c_arc3D));
poms_find_resk_lonb(pt.alpha,&CA[0]);
poms_find_resk_lonb(pt.beta ,&CA[1]);
poms_find_resk_lonb(pt.gamma,&CA[2]);
poms_find_resk_lonb(pt.delta,&CA[3]);
fewg_find_duvk(probe,S1,S2,&omega);
culm_unit_peks(S1.zent,S2.zent,&nrm);

diff=(double *)malloc(4*sizeof(double));
sml=LARGE_NUMBER;
for(i=0;i<4;i++)
	{sp=rocv_scal_toqc(nrm,CA[i].nrml);
	fb=fabs(sp);
	diff[i]=fabs(fb-1.0);
	if(diff[i]<sml)
		{sml=diff[i];
		q1=i;
		}
	}
poms_find_resk_lonb(CA[q1],&C1);

sml=LARGE_NUMBER;
for(i=0;i<4;i++)if(i!=q1)
	{if(diff[i]<sml)
		{sml=diff[i];
		q2=i;
		}
	}
poms_find_resk_lonb(CA[q2],&C2);

p1=-1;  p2=-1;
for(i=0;i<4;i++)
if((i!=q1)&&(i!=q2))
	{if(p1==-1)		cs=1;
	else			cs=2;
	if(cs==1)
		{poms_find_resk_lonb(CA[i],&D_1);
		p1=i;
		}
	if(cs==2)
		{poms_find_resk_lonb(CA[i],&D_2);
		break;
		}
	}

nb=zolf_find_fudp_miwf(N,C1,C2,D_1,D_2,omega,PT,max_PT);
free(diff);
free(CA);
return nb;
}


void zifb_unif_tolp(c_arc3D CA,int N,
c_arc3D *C,point *sep)
{int n,i;
double a,b,t,step;
n=N-1;
qirp_inte_ligr(CA,&a,&b);
step=(b-a)/(double)N;
for(i=1;i<=N-1;i++)
	{t=a+(double)i*step;
	nehl_eval_segt(CA,t,&sep[i-1]);
	}
sahf_spli_gehq(CA,sep,n,C);
}


double milc_dist_leqk(c_arc3D C1,
c_arc3D C2)
{int dummy;
double d_s,d_t,dis1,dis2,dis,res,diff_md;
point mid1,mid2;
d_s=cijv_dist_laph(C1.begn,C2,&dummy);
d_t=cijv_dist_laph(C1.term,C2,&dummy);
if(d_s<d_t)	dis1=d_t;
else		dis1=d_s;

d_s=cijv_dist_laph(C2.begn,C1,&dummy);
d_t=cijv_dist_laph(C2.term,C1,&dummy);
if(d_s<d_t)	dis2=d_t;
else		dis2=d_s;
if(dis1>dis2)	dis=dis1;
else			dis=dis2;

renw_midp_mocw(C1,&mid1);
renw_midp_mocw(C2,&mid2);
diff_md=wodt_dist_gilq(mid1,mid2);
res=dis+diff_md;
return res;
}


double pesc_dist_tejw(pt_tor PT,
c_arc3D C)
{int i;
double dis,res;
c_arc3D *CA;
CA=(c_arc3D *)malloc(4*sizeof(c_arc3D));
poms_find_resk_lonb(PT.alpha,&CA[0]);
poms_find_resk_lonb(PT.beta ,&CA[1]);
poms_find_resk_lonb(PT.gamma,&CA[2]);
poms_find_resk_lonb(PT.delta,&CA[3]);
res=LARGE_NUMBER;
for(i=0;i<4;i++)
	{dis=milc_dist_leqk(CA[i],C);
	if(dis<res)
		res=dis;
	}
free(CA);
return res;
}



int wanz_inci_vufq(set_arcs SA,pt_tor PT,
c_arc3D *C)
{int i,N,q,id;
double err,sml;
N=SA.ar_grs;
if(N<1)
	{fprintf(tmpout,"Unable to find arc\n");
	exit(0);
	}
sml=LARGE_NUMBER;
for(i=0;i<N;i++)
	{err=pesc_dist_tejw(PT,SA.C[i]);
	if(err<sml)
		{sml=err;
		q=i;
		}
	}
poms_find_resk_lonb(SA.C[q],C);
id=q;
return id;
}


double vogp_find_wafk_piws(trmsrf surf,c_curve cc,int comp,
c_arc3D C)
{double a,b,m,res;
point MD_C,MD_N,temp;
kehf_inte_recn(cc,comp,&a,&b);
m=0.5*(a+b);
novc_eval_vokn(cc,m,&temp);
wolf_eval_murg(surf,temp.absi,temp.ordo,&MD_N);
renw_midp_mocw(C,&MD_C);
res=wodt_dist_gilq(MD_C,MD_N);
return res;
}



double jicw_find_dofg_nugw(trmsrf ts,c_curve cc,
c_arc3D CA,int *idx)
{int nb_ct,i,nx;
double res,dis,sml,D;
point *sep,A,B;
nb_ct=cc.N;
sep=(point *)malloc(nb_ct*sizeof(point));
sofl_segm_salc(ts,cc,sep);
sml=LARGE_NUMBER;
for(i=0;i<nb_ct;i++)
	{nx=i+1;
	if(nx==nb_ct)
		nx=0;
	getf_find_rogc_todj(sep[i],&A);
	getf_find_rogc_todj(sep[nx],&B);
	D=vogp_find_wafk_piws(ts,cc,i,CA);
	dis=cutj_dist_rulb(A,B,CA);
	if(D<sml)
		{sml=D;
		*idx=i;
		res=dis;
		}
	}
free(sep);
return res;
}



double cewl_find_ziwb_hutn(trmsrf ts,c_curve cc,
c_arc3D CA,int *idx)
{int nb_ct,i,nx;
double res,dis;
point *sep,A,B;
nb_ct=cc.N;
sep=(point *)malloc(nb_ct*sizeof(point));
sofl_segm_salc(ts,cc,sep);
res=LARGE_NUMBER;
for(i=0;i<nb_ct;i++)
	{nx=i+1;
	if(nx==nb_ct)
		nx=0;
	getf_find_rogc_todj(sep[i],&A);
	getf_find_rogc_todj(sep[nx],&B);
	dis=cutj_dist_rulb(A,B,CA);
	if(dis<res)
		{res=dis;
		*idx=i;
		}
	}
free(sep);
return res;
}


double figh_find_dotq_buvr(trmsrf ts,c_curve cc,
c_arc3D CA,int *idx)
{int nb_ct;
double res;
nb_ct=cc.N;
if(nb_ct==2)
	res=jicw_find_dofg_nugw(ts,cc,CA,idx);
else
	res=cewl_find_ziwb_hutn(ts,cc,CA,idx);
return res;
}



double qevd_find_guwk_muck(trmsrf ts,c_arc3D CA,int *bound,int *idx)
{int nin,id,i;
double err,res;

res=figh_find_dotq_buvr(ts,ts.cc,CA,&id);
*idx=id;
*bound=-1;

nin=ts.nb_inner;
for(i=0;i<nin;i++)
	{err=figh_find_dotq_buvr(ts,ts.inner[i],CA,&id);
	if(err<res)
		{*bound=i;
		*idx=id;
		res=err;
		}
	}
return res;
}



void qirw_comp_nuwd(c_arc3D CA,int p,trmsrf *surf1,
 int nb_surf1,int *supp,int *surf_id,int *bound_id,int *comp_id)
{int nb_max=10,i,k,*upd,n_upd,w;
int bound,idx;
double sml,dis;
upd=(int *)malloc(nb_max*sizeof(int));
k=0;
for(i=0;i<nb_surf1;i++)
if(supp[i]==p)
	{upd[k]=i;
	k++;
	if(k>=nb_max)
		{fprintf(tmpout,"nb_max is had\n");
		exit(0);
		}
	}
n_upd=k;

sml=LARGE_NUMBER;
for(i=0;i<n_upd;i++)
	{w=upd[i];
	dis=qevd_find_guwk_muck(surf1[w],CA,&bound,&idx);
	if(dis<sml)
		{sml=dis;
		*surf_id=w;
		*bound_id=bound;
		*comp_id=idx;
		}
	}
free(upd);
}


double kenh_find_zocq_jegm(ns_curv nc,point beg,
int *ort,point *new_beg)
{int n;
double dis_st,dis_tr,res;
parm ST,TR,beg2D;
n=nc.n;
ST.u=nc.d[0].absi;
ST.v=nc.d[0].ordo;
TR.u=nc.d[n].absi;
TR.v=nc.d[n].ordo;
beg2D.u=beg.absi;
beg2D.v=beg.ordo;
dis_st=pufv_dist_mekq(beg2D,ST);
dis_tr=pufv_dist_mekq(beg2D,TR);
if(dis_st<dis_tr)
	{res=dis_st;
	getf_find_rogc_todj(nc.d[n],new_beg);
	*ort=+1;
	}
else
	{res=dis_tr;
	getf_find_rogc_todj(nc.d[0],new_beg);
	*ort=-1;
	}
return res;
}



void lemz_upda_jadh(trmsrf srf,c_arc3D *C,
int nb_segm,c_curve *cc,int id)
{int i,j,k,nb_cp,new_cp,ort;
int *seq,*exc,q,*or;
double a,b,dis,sml;
point beg,new_beg,cur_beg;
ns_curv *nc,*temp;
prop_n_curv pnc;

if(srf.type!=3)
	{fprintf(tmpout,"trimmed surface must have type=3 here\n");
	exit(0);
	}
pnc.n=4;
pnc.k=3;
nc=(ns_curv *)malloc(nb_segm*sizeof(ns_curv));
for(i=0;i<nb_segm;i++)
	{foks_allo_vukp(pnc,&nc[i]);
	dolj_curv_kacq(srf.ts,C[i],&nc[i]);
	}

nb_cp=cc->N;
if(nb_cp-1+nb_segm>=MAXCOMP)
	{fprintf(tmpout,"1.  MAXCOMP=%d is reached\n",MAXCOMP);
	exit(0);
	}
kehf_inte_recn(*cc,id,&a,&b);
novc_eval_vokn(*cc,a,&beg);
exc=(int *)malloc(nb_segm*sizeof(int));
for(i=0;i<nb_segm;i++)
	exc[i]=0;
seq=(int *)malloc(nb_segm*sizeof(int));
or=(int *)malloc(nb_segm*sizeof(int));
for(j=0;j<nb_segm;j++)
	{sml=LARGE_NUMBER;
	for(i=0;i<nb_segm;i++)if(exc[i]==0)
		{dis=kenh_find_zocq_jegm(nc[i],beg,&ort,&new_beg);
		if(dis<sml)
			{sml=dis;
			q=i;
			seq[j]=i;
			or[j]=ort;
			getf_find_rogc_todj(new_beg,&cur_beg);
			}
		}
	getf_find_rogc_todj(cur_beg,&beg);
	exc[q]=+1;
	}
free(exc);

temp=(ns_curv *)malloc(nb_segm*sizeof(ns_curv));
for(i=0;i<nb_segm;i++)
	foks_allo_vukp(pnc,&temp[i]);
for(i=0;i<nb_segm;i++)
	{q=seq[i];
	if(or[i]==+1)
		zobm_find_wumq_kihf(nc[q],&temp[i]);
	if(or[i]==-1)
		colw_inve_pelj(nc[q],&temp[i]);
	}
free(seq);
free(or);
for(i=0;i<nb_segm;i++)
	zobm_find_wumq_kihf(temp[i],&nc[i]);
for(i=0;i<nb_segm;i++)
	newt_dest_lefq(pnc,&temp[i]);
free(temp);

new_cp=nb_cp-1+nb_segm;
temp=(ns_curv *)malloc(new_cp*sizeof(ns_curv));
for(i=0;i<new_cp;i++)
	foks_allo_vukp(pnc,&temp[i]);
k=0;
for(i=0;i<id;i++)
	{zobm_find_wumq_kihf(cc->nc[i],&temp[k]);
	k++;
	}
for(i=0;i<nb_segm;i++)
	{zobm_find_wumq_kihf(nc[i],&temp[k]);
	k++;
	}
for(i=id+1;i<nb_cp;i++)
	{zobm_find_wumq_kihf(cc->nc[i],&temp[k]);
	k++;
	}
for(i=0;i<nb_segm;i++)
	newt_dest_lefq(pnc,&nc[i]);
free(nc);

for(i=0;i<new_cp;i++)
	{zobm_find_wumq_kihf(temp[i],&cc->nc[i]);
	cc->type[i]=2;
	}
cc->N=new_cp;
cc->nnc=new_cp;
cc->nca=0;
cc->nle=0;
for(i=0;i<new_cp;i++)
	newt_dest_lefq(pnc,&temp[i]);
free(temp);
}



void cibj_upda_vutd(int p,c_arc3D CA,c_arc3D *C,int nb_segm,
trmsrf *surf1,int nb_surf1,int *supp)
{int s_id,bd_id,cp_id;
qirw_comp_nuwd(CA,p,surf1,nb_surf1,supp,&s_id,&bd_id,&cp_id);
if(bd_id==-1)
	lemz_upda_jadh(surf1[s_id],C,nb_segm,&surf1[s_id].cc,cp_id);
else if(bd_id>=0)
	lemz_upda_jadh(surf1[s_id],C,nb_segm,&surf1[s_id].inner[bd_id],cp_id);

}



void kemp_upda_gerc(atom *A,int nb_segm,trmsrf *surf2,
blend_cpx BC,int z,set_arcs *SA,trmsrf *surf1,
int nb_surf1,int *supp)
{int *p,i,j,id,w,nb,n;
point *sep;
c_arc3D CA,*C;
w=BC.BT[z].trim_idx;
p=(int *)malloc(3*sizeof(int));
C=(c_arc3D *)malloc(nb_segm*sizeof(c_arc3D));
p[1]=BC.BT[z].sph_idx1;
p[2]=BC.BT[z].sph_idx2;
n=nb_segm-1;
sep=(point *)malloc(n*sizeof(point)); 
for(i=1;i<=2;i++)
	{id=wanz_inci_vufq(SA[p[i]],surf2[w].pt,&CA);
	zifb_unif_tolp(CA,nb_segm,C,sep);
	nb=SA[p[i]].ar_grs;
	poms_find_resk_lonb(C[0],&SA[p[i]].C[id]);
	for(j=1;j<nb_segm;j++)
		{if(nb+j-1>=MAX_ARCS)
			{fprintf(tmpout,"MAX_ARCS is reached\n");
			exit(0);
			}
		poms_find_resk_lonb(C[j],&SA[p[i]].C[nb+j-1]);
		}
	SA[p[i]].ar_grs=nb+nb_segm-1;
		
	cibj_upda_vutd(p[i],CA,C,nb_segm,surf1,nb_surf1,supp);
	}
free(sep);
free(C);
free(p);
}


int hefk_apen_sikf(double probe,int N,sphere *S,int nb_sph,
trmsrf *surf1,int nb_surf1,int *supp,trmsrf *surf2,
int nb_surf2,blend_cpx *BC,set_arcs *SA,int max_surf2)
{int i,j,w,p1,p2,k,newnb;
int nb_blend,nb,k_bl,max_PT;
pt_tor *PT;
nb_blend=BC->bt_grs;
if(nb_blend>=MAX_BLEND_TOR)
	{fprintf(tmpout,"MAX_BLEND_TOR has been exceeded\n");
	exit(0);
	}
max_PT=N+4;
PT=(pt_tor *)malloc(max_PT*sizeof(pt_tor));
k=nb_surf2;	k_bl=nb_blend;
for(i=0;i<nb_blend;i++)
	{if(verbose_variable==VERBOSE)
		{if((i% 20==0)||(i==nb_blend-1))
			fprintf(tmpout,"Blend=%d / %d  cur loc=%d\n",i,nb_blend-1,k);
		}
	w =BC->BT[i].trim_idx;
	p1=BC->BT[i].sph_idx1;
	p2=BC->BT[i].sph_idx2;
	if((w>=nb_surf2)||(w<0))
		{fprintf(tmpout,"Access beyond existing surf2:  w=%d\n",w);
		exit(0);
		}
	nb=setw_deco_wipz(S[p1],S[p2],probe,N,surf2[w].pt,PT,max_PT);
	if(nb>=2)
		{kemp_upda_gerc(S,nb,surf2,*BC,i,SA,surf1,nb_surf1,supp);
		
		romh_find_cont_qucr(PT[0],&surf2[w].pt);
		
		BC->BT[i].trim_idx=w;
		for(j=1;j<nb;j++)
			{if(k>=max_surf2)
				{fprintf(tmpout,"current max=%d\n",max_surf2);
				fprintf(tmpout,"Max nb trimmed surfaces is reached\n");
				exit(0);
				}
			surf2[k].type=4;
			surf2[k].boundary=1;
			romh_find_cont_qucr(PT[j],&surf2[k].pt);
			sart_rect_jamc(0.0,1.0,0.0,1.0,&surf2[k].cc);
			surf2[k].nb_inner=0;
			if(k_bl>=MAX_BLEND_TOR)
				{fprintf(tmpout,"MAX_BLEND_TOR is reached\n");
				exit(0);
				}
			BC->BT[k_bl].trim_idx=k;
			BC->BT[k_bl].sph_idx1=p1;
			BC->BT[k_bl].sph_idx2=p2;
			k++;	k_bl++;
			if(k>=max_surf2)
				{fprintf(tmpout,"current max=%d\n",max_surf2);
				fprintf(tmpout,"Max nb trimmed surfaces is reached\n");
				exit(0);
				}
			}
		}
	}
newnb=k;
free(PT);
BC->bt_grs=k_bl;
kejh_fill_zogq(S,surf2,nb_sph,BC);
return newnb;
}


