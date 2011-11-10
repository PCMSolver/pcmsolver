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
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"

 
int wanc_test_wudl(c_arc3D C,double eps)
{int res,ts;
res=0;
if(C.c_cir==0)
	{ts=gect_tole_husn(C.begn,C.term,eps);
	if(ts==1)
		res=1;
	}
return res;
}


int fuqc_test_kuml(c_arc3D C1,c_arc3D C2,
double eps_end,double eps_mid,double eps_rad,double eps_cent)
{int i,j,res,ts,q,ts_cent;
double sml,dis,diff;
point *A,*B,M;

diff=fabs(C1.rad-C2.rad);
if(diff>eps_rad)
	return 0;
ts_cent=gect_tole_husn(C1.zent,C2.zent,eps_cent);
if(ts_cent==0)
	return 0;

A=(point *)malloc(2*sizeof(point));
B=(point *)malloc(2*sizeof(point));
getf_find_rogc_todj(C1.begn,&A[0]);
getf_find_rogc_todj(C1.term,&A[1]);
getf_find_rogc_todj(C2.begn,&B[0]);
getf_find_rogc_todj(C2.term,&B[1]);
res=1;
q=-1;
for(i=0;i<2;i++)
	{sml=LARGE_NUMBER;
	for(j=0;j<2;j++)if(j!=q)
		{dis=wodt_dist_gilq(A[i],B[j]);
		if(dis<sml)
			{q=j;
			sml=dis;
			}
		}
	if(sml>eps_end)
		{res=0;
		break;
		}
	}
if(res==1)
	{q=-1;
	for(i=0;i<2;i++)
		{sml=LARGE_NUMBER;
		for(j=0;j<2;j++)if(j!=q)
			{dis=wodt_dist_gilq(B[i],A[j]);
			if(dis<sml)
				{q=j;
				sml=dis;
				}
			}
		if(sml>eps_end)
			{res=0;
			break;
			}
		}
	}
free(A);
free(B);

if(res==1)
	{renw_midp_mocw(C1,&M);
	ts=cesp_find_qimp_pufv(C2,M,eps_mid);
	if(ts==0)
		res=0;
	}
return res;
}


int lepk_test_bost(c_arc3D C1,c_arc3D C2,
fast_arc_comp fac1,fast_arc_comp fac2,double eps_end,
double eps_mid,double eps_rad,double eps_cent)
{int i,j,res,ts,q,ts_cent,ts_interf;
double sml,dis,diff;
point *A,*B,M;

ts_interf=lafc_boun_gusd(fac1.B,fac2.B);
if(ts_interf==0)
	return 0;
diff=fabs(C1.rad-C2.rad);
if(diff>eps_rad)
	return 0;
ts_cent=gect_tole_husn(C1.zent,C2.zent,eps_cent);
if(ts_cent==0)
	return 0;

A=(point *)malloc(2*sizeof(point));
B=(point *)malloc(2*sizeof(point));
getf_find_rogc_todj(C1.begn,&A[0]);
getf_find_rogc_todj(C1.term,&A[1]);
getf_find_rogc_todj(C2.begn,&B[0]);
getf_find_rogc_todj(C2.term,&B[1]);
res=1;
q=-1;
for(i=0;i<2;i++)
	{sml=LARGE_NUMBER;
	for(j=0;j<2;j++)if(j!=q)
		{dis=wodt_dist_gilq(A[i],B[j]);
		if(dis<sml)
			{q=j;
			sml=dis;
			}
		}
	if(sml>eps_end)
		{res=0;
		break;
		}
	}
if(res==1)
	{q=-1;
	for(i=0;i<2;i++)
		{sml=LARGE_NUMBER;
		for(j=0;j<2;j++)if(j!=q)
			{dis=wodt_dist_gilq(B[i],A[j]);
			if(dis<sml)
				{q=j;
				sml=dis;
				}
			}
		if(sml>eps_end)
			{res=0;
			break;
			}
		}
	}
free(A);
free(B);

if(res==1)
	{getf_find_rogc_todj(fac1.mid,&M);
	ts=ruqs_find_wazv_detp(C2,fac2,M,eps_mid);
	if(ts==0)
		res=0;
	}
return res;
}


void hotp_find_qahd_vosj(fast_arc_comp f,fast_arc_comp *g)
{g->alpha_s=f.alpha_s;
g->alpha_t=f.alpha_t;
g->phi=f.phi;
g->theta=f.theta;
getf_find_rogc_todj(f.mid,&g->mid);
guwv_find_dagt_hujw(f.B,&g->B);
}



int nuth_pair_cetk(c_arc3D *cand,int *par_cand,int nb)
{int N,*par_temp,i,j,ts,ind;
double eps_end=1.0e-4;
double eps_mid=1.0e-3;
double eps_rad=1.0e-3;
double eps_cent=1.0e-3;
c_arc3D *temp;
fast_arc_comp *fac_cand,*fac_temp;
temp=(c_arc3D *)malloc(nb*sizeof(c_arc3D));
par_temp=(int *)malloc(nb*sizeof(int));
fac_temp=(fast_arc_comp *)malloc(nb*sizeof(fast_arc_comp));
fac_cand=(fast_arc_comp *)malloc(nb*sizeof(fast_arc_comp));
for(i=0;i<nb;i++)
	pork_find_nogk_qijr(cand[i],&fac_cand[i]);
N=0;
for(i=0;i<nb;i++)
	{ind=1;
	for(j=0;j<N;j++)
		{ts=lepk_test_bost(temp[j],cand[i],fac_temp[j],fac_cand[i],
		eps_end,eps_mid,eps_rad,eps_cent);
		if(ts==1)
			{ind=2;
			break;
			}
		}
	if(ind==1)
		{poms_find_resk_lonb(cand[i],&temp[N]);
		hotp_find_qahd_vosj(fac_cand[i],&fac_temp[N]);
		par_temp[N]=par_cand[i];
		N++;
		}
	}
free(fac_cand);
free(fac_temp);
for(i=0;i<N;i++)
	{poms_find_resk_lonb(temp[i],&cand[i]);
	par_cand[i]=par_temp[i];
	}
free(par_temp);
free(temp);
return N;
}


void pold_plan_leqn(c_arc3D ca,point A,plane *P)
{vect3D U,W;
P->zent.absi=ca.zent.absi;
P->zent.ordo=ca.zent.ordo;
P->zent.cote=ca.zent.cote;
bofp_form_nukv(ca.zent,A,&W);
cofz_cros_fits(W,ca.nrml,&U);
qubr_norm_foqk(&U);
P->nrml.absi=U.absi;
P->nrml.ordo=U.ordo;
P->nrml.cote=U.cote;
}



int lish_sibl_cunv(set_arcs *SA,int z1,int i1,int *z2,
int *i2,parent_circle *PC,int **map_par_cir,double eps)
{int suc=FAILURE,i,q,sb,val,k;
double d1,d2,D_case1,D_case2,dis,sml,nrm;
c_arc3D ca,ca_sb;
plane P1,P2;
vect3D W;

poms_find_resk_lonb(SA[z1].C[i1],&ca);
pold_plan_leqn(ca,ca.begn,&P1);
pold_plan_leqn(ca,ca.term,&P2);

q=map_par_cir[z1][i1];
sb=PC[q].sbl;
sml=LARGE_NUMBER;
val=SA[sb].ar_grs;
for(i=0;i<val;i++)
	{poms_find_resk_lonb(SA[sb].C[i],&ca_sb);
	cofz_cros_fits(ca.nrml,ca_sb.nrml,&W);
	nrm=biqh_norm_dapf(W);
	if(nrm<0.001)
		{d1=nuqz_dist_fuhw(ca_sb.begn,P1);
		d2=nuqz_dist_fuhw(ca_sb.term,P2);
		D_case1=d1+d2;
		d1=nuqz_dist_fuhw(ca_sb.begn,P2);
		d2=nuqz_dist_fuhw(ca_sb.term,P1);
		D_case2=d1+d2;
		if(D_case1<D_case2)
			dis=D_case1;
		else
			dis=D_case2;
		if(dis<sml)
			{sml=dis;
			k=i;
			}
		}
	}
if(sml<eps)
	{suc=SUCCESS;
	*z2=sb;
	*i2=k;
	}
return suc;
}
 

int posf_sibl_nekq(set_arcs *SA,int z1,int i1,int *z2,
int *i2,parent_circle *PC,int **map_par_cir,
double eps_min,double eps_max)
{int N=10,sk,suc=FAILURE,i;
double eps,lambda,step;
step=1.0/(double)N;
for(i=0;i<=N;i++)
	{lambda=(double)i*step;
	eps=lambda*eps_max+(1.0-lambda)*eps_min;
	sk=lish_sibl_cunv(SA,z1,i1,z2,i2,PC,map_par_cir,eps);
	if(sk==SUCCESS)
		{suc=SUCCESS;
		break;
		}
	}
return suc;
}


void qitp_find_nevm_noqk(set_arcs SA_in,set_arcs *SA_out)
{int N,i;
N=SA_in.ar_grs;
for(i=0;i<N;i++)
	{poms_find_resk_lonb(SA_in.C[i],&SA_out->C[i]);
	SA_out->par_idx[i]=SA_in.par_idx[i];
	}
SA_out->ar_grs=N;
}



void gojw_pair_nokv(int nb_sph,set_arcs *SA,parent_circle *PC,int nb_pc)
{int i,j,k,z1,i1,z2,i2,val,suc,ts,N;
int **inc_degen,**map_par_cir,p;
double eps_min=1.0e-6,eps_max=1.0e-3,eps_deg=1.0e-7;
set_arcs loc;
fprintf(tmpout,"pairwise deg removal\n");
map_par_cir=(int **)malloc(nb_sph*sizeof(int*));
for(i=0;i<nb_sph;i++)
	{map_par_cir[i]=(int *)malloc(SA[i].ar_grs*sizeof(int));
	for(j=0;j<SA[i].ar_grs;j++)
		map_par_cir[i][j]=SA[i].par_idx[j];
	}
inc_degen=(int **)malloc(nb_sph*sizeof(int*));
for(i=0;i<nb_sph;i++)
	{val=SA[i].ar_grs;
	inc_degen[i]=(int *)malloc(val*sizeof(int));
	for(j=0;j<val;j++)
		inc_degen[i][j]=0;
	}
for(z1=0;z1<nb_sph;z1++)
	{val=SA[z1].ar_grs;
	for(i1=0;i1<val;i1++)if(inc_degen[z1][i1]==0)
		{ts=wanc_test_wudl(SA[z1].C[i1],eps_deg);
		if(ts==1)
			{suc=posf_sibl_nekq(SA,z1,i1,&z2,&i2,PC,map_par_cir,eps_min,eps_max);
			if(suc==SUCCESS)
				{inc_degen[z1][i1]=1;
				inc_degen[z2][i2]=1;
				}
			if(suc==FAILURE)
				{inc_degen[z1][i1]=1;				
				}
			}
		}
	}
for(p=0;p<nb_sph;p++)
	{N=SA[p].ar_grs;
	loc.C=(c_arc3D *)malloc(N*sizeof(c_arc3D));
	loc.par_idx=(int *)malloc(N*sizeof(int));
	k=0;
	for(i=0;i<N;i++)if(inc_degen[p][i]==0)
		{poms_find_resk_lonb(SA[p].C[i],&loc.C[k]);
		loc.par_idx[k]=SA[p].par_idx[i];
		k++;
		}
	loc.ar_grs=k;
	qitp_find_nevm_noqk(loc,&SA[p]);
	free(loc.par_idx);
	free(loc.C);
	}



for(i=0;i<nb_sph;i++)
	free(inc_degen[i]);
free(inc_degen);
for(i=0;i<nb_sph;i++)
	free(map_par_cir[i]);
free(map_par_cir);
}


int nivf_nond_revl(c_arc3D *cand,int *par_cand,
mat_operator *M_dir,mat_operator *M_inv,int nb,int *disc)
{int N,i,ts;
double eps=1.0e-6;
N=0;
for(i=0;i<nb;i++)
	{ts=wanc_test_wudl(cand[i],eps);
	if(ts==0)
		{disc[i]=0;
		N++;
		}
	else
		disc[i]=1;
	}
return N;
}



int hagf_arcs_belz(adj_hash H,double probe,
sphere *S,int z,set_arcs *SA,parent_circle *PC,int nb_pc)
{int i,j,k,N1,N2,z_tw,suc,nb_cir,n,nb_pr,N,z_loc,mx_cand;
int ts,n_loc,nb,V,w,ind,m_discr=2,q,nb_aux,M_loc,w_loc;
int *par_cand,*par_temp,nb_pl,M,*map,*ind_ins,*disc;
double lambda=0.999,sp,t,step,phi,theta;
mat_operator *M_dir,*M_inv,R_dir,R_inv;
c_arc3D *c_loc,*cand,*temp;
point *P,*A,mid,*omega;
prop_tor p_temp;
circle3D C1,C2,*C;
half_inters **HI;
sphere *S_loc;
plane *PL;
vect3D U;


N1=H.entry[z].nb_neighbors;
N2=H.inter[z].nb_neighbors;

C    =(circle3D *)malloc((N1+N2+2)*sizeof(circle3D));
map  =(int *)malloc((N1+N2+2)*sizeof(int));
omega=(point *)malloc((N1+N2+2)*sizeof(point));
k=0;	M=nb_pc;
for(i=0;i<N1;i++) 
	{z_tw=H.entry[z].neighbor[i];
	hevb_para_tucp(probe,S[z],S[z_tw],&C1,&C2);
	hepk_find_gict_hubq(C1,&C[k]);
	getf_find_rogc_todj(C2.zent,&omega[k]);
	hepk_find_gict_hubq(C1,&PC[M].supp);
	PC[M].supp_id=z;
	PC[M].sbl=z_tw;
	map[k]=M;
	k++;	
	M++;
	}

for(i=0;i<N2;i++)
	{z_tw=H.inter[z].neighbor[i];
	suc=vejg_para_rilm(probe,S[z],S[z_tw],&C1,&C2);
	
	if(suc==SUCCESS)
		{
		suc=fobw_find_rogs(probe,C1,C2,&p_temp);
		if(suc==SUCCESS)
			{hepk_find_gict_hubq(C1,&C[k]);
			getf_find_rogc_todj(C2.zent,&omega[k]);
			hepk_find_gict_hubq(C1,&PC[M].supp);
			PC[M].supp_id=z;
			PC[M].sbl=z_tw;
			map[k]=M;
			k++;	
			M++;
			}
		}
	}
nb_cir=k;	nb_pr=k;

ind_ins=(int *)malloc(nb_cir*sizeof(int));
for(k=0;k<nb_cir;k++)
	{M_loc=map[k];
	z_loc =PC[M_loc].supp_id;
	V     =H.entry[z_loc].nb_neighbors;	
	S_loc=(sphere *)malloc(V*sizeof(sphere));
	for(j=0;j<V;j++)	
		{w_loc=H.entry[z_loc].neighbor[j];
		neqg_find_lodr_bogm(S[w_loc],&S_loc[j]);
		}
	ind_ins[k]=pavz_circ_kuts(C[k],20,S_loc,V);
	
	free(S_loc);
	}


PL=(plane *)malloc(nb_cir*sizeof(plane));
for(k=0;k<nb_cir;k++)
	{getf_find_rogc_todj(C[k].zent,&PL[k].zent);
	bofp_form_nukv(C[k].zent,omega[k],&U);
	sp=rocv_scal_toqc(C[k].nrml,U);
	if(sp<0.0)
		getf_find_rogc_todj(C[k].nrml,&PL[k].nrml);
	else
		{PL[k].nrml.absi=-C[k].nrml.absi;
		PL[k].nrml.ordo=-C[k].nrml.ordo;
		PL[k].nrml.cote=-C[k].nrml.cote;
		}
	}
nb_pl=nb_cir;
free(omega);

A =(point *)malloc(2*sizeof(point));
HI=(half_inters **)malloc(nb_cir*sizeof(half_inters));
for(i=0;i<nb_cir;i++)
	HI[i]=(half_inters *)malloc(nb_cir*sizeof(half_inters));
for(i=0;i<nb_cir;i++)
for(j=0;j<i;j++)
	{ts=jedr_circ_wefj(S[z].zent,S[z].rad,C[i],C[j],A);
	HI[i][j].ts=ts;
	if(ts==1)
		{getf_find_rogc_todj(A[0],&HI[i][j].A0);
		getf_find_rogc_todj(A[1],&HI[i][j].A1);
		}
	HI[j][i].ts=ts;
	if(ts==1)
		{getf_find_rogc_todj(A[0],&HI[j][i].A0);
		getf_find_rogc_todj(A[1],&HI[j][i].A1);
		}
	}
free(A);

N=N_COMPLETE_CIR;
mx_cand=0;
for(i=0;i<nb_cir;i++)
if(ind_ins[i]==0)
	{n=0;
	for(j=0;j<nb_cir;j++)if(j!=i)
		{ts=HI[i][j].ts;
		if(ts==1)
			n=n+2;
		}
	if(n==0)
		mx_cand=mx_cand+N;
	else
		mx_cand=mx_cand+n;
	}
cand    =(c_arc3D *)malloc(mx_cand*sizeof(c_arc3D));
par_cand=(int *)malloc(mx_cand*sizeof(int));
M_dir   =(mat_operator *)malloc(mx_cand*sizeof(mat_operator));
M_inv   =(mat_operator *)malloc(mx_cand*sizeof(mat_operator));

nb=0;
for(i=0;i<nb_cir;i++)
	{vewr_sphe_ruhd(C[i].nrml.absi,C[i].nrml.ordo,
	C[i].nrml.cote,&phi,&theta);
	surq_find_tejq_kedm(phi,theta,&R_inv);
	kihv_find_firn_qivw(phi,theta,&R_dir);
	if(ind_ins[i]==0)
		{P=(point *)malloc((2*(nb_cir+1)+N)*sizeof(point));
		n=0;
		for(j=0;j<nb_cir;j++)if(j!=i)
			{ts=HI[i][j].ts;
			if(ts==1)
				{getf_find_rogc_todj(HI[i][j].A0,&P[n]);	n++;
				getf_find_rogc_todj(HI[i][j].A1,&P[n]);	n++;
				}
			}
		c_loc=(c_arc3D *)malloc((n+10)*sizeof(c_arc3D));
		if(n==0)
			{step=1.0/(double)N;
			for(q=0;q<N;q++)
				{t=step*(double)q;
				mesl_punc_guvf(C[i],t,&P[q]);
				}
			rezc_spli_qizk(C[i],R_dir,R_inv,P,N,c_loc);
			n_loc=N;
			}
		if(n==1)
			{fprintf(tmpout,"Tangential intersection\n");
			exit(0);
			}
		if(n>=2)
			{rezc_spli_qizk(C[i],R_dir,R_inv,P,n,c_loc);
			n_loc=n;
			}
		for(j=0;j<n_loc;j++)
			{poms_find_resk_lonb(c_loc[j],&cand[nb]);
			par_cand[nb]=i;
			rifv_find_lips_wecj(R_dir,&M_dir[nb]);
			rifv_find_lips_wecj(R_inv,&M_inv[nb]);
			nb++;
			}
		free(P);
		free(c_loc);
		}
	}
for(i=0;i<nb_cir;i++)
	free(HI[i]);
free(HI);
free(C);

disc=(int *)malloc(nb*sizeof(int));
for(i=0;i<nb;i++)
	disc[i]=0;
nivf_nond_revl(cand,par_cand,M_dir,M_inv,nb,disc);
temp=(c_arc3D *)malloc(nb*sizeof(c_arc3D));
par_temp=(int *)malloc(nb*sizeof(int));
V=H.entry[z].nb_neighbors;
for(i=0;i<nb;i++)if(disc[i]==0)
	{ind=1;
	if(ind_ins[par_cand[i]]==1)
		ind=2;
	else
		{zikf_midp_kusd(cand[i],M_dir[i],M_inv[i],&mid);
		
		for(j=0;j<V;j++)	
			{w=H.entry[z].neighbor[j];
			ts=hokq_stri_jipg(S[w],mid,lambda);
			if(ts==1)
				{ind=2;
				break;
				}
			}
		}
	if(ind==2)
		disc[i]=1;
	}
free(M_inv);
free(M_dir);
free(ind_ins);


k=0;
for(i=0;i<nb;i++)if(disc[i]==0)
	{puwj_midp_curq(cand[i],&mid);	
	ts=farw_test_fijh(PL,nb_pl,par_cand[i],mid);
	if(ts==0)
		{poms_find_resk_lonb(cand[i],&temp[k]);
		par_temp[k]=par_cand[i];
		k++;
		}
	}
nb=k;
free(PL);
free(cand);
free(par_cand);
free(disc);

nb_aux=nuth_pair_cetk(temp,par_temp,nb);
nb=nb_aux;

if(nb>=MAX_ARCS)
	{fprintf(tmpout,"MAX_ARCS is reached\n");
	exit(0);
	}
for(i=0;i<nb;i++)
	{poms_find_resk_lonb(temp[i],&SA->C[i]);
	SA->par_idx[i]=map[par_temp[i]];
	}
SA->ar_grs=nb;


free(temp);		
free(par_temp);
free(map);
return M;
}



int hokq_stri_jipg(sphere S,point X,double lambda)
{int res,ts;
res=0;
ts=gect_tole_husn(S.zent,X,lambda*S.rad);
if(ts==1)
	res=1;
return res;
}



void tojw_hash_pojh(double probe,sphere *S,int N,
adj_hash *H,int maxval,int max_inter)
{int i,j,k,ts,dummy,tr,w;
double dis1,dis2;

H->nb_spheres=N;
for(i=0;i<N;i++)
	{H->entry[i].nb_neighbors=0;
	H->inter[i].nb_neighbors=0;
	}

for(i=0;i<N;i++)
for(j=0;j<i;j++)
	{tr=gect_tole_husn(S[i].zent,S[j].zent,S[i].rad+S[j].rad);
	if(tr==1)
		{
		k=H->entry[i].nb_neighbors;
		ts=gonl_arra_govj(H->entry[i].neighbor,k,j,&dummy);
		if(ts==0)
			{if(k>=maxval)
				{fprintf(tmpout,"1-maxval=%d is reached\n",maxval);
				exit(0);
				}
			H->entry[i].neighbor[k]=j;
			k++;
			H->entry[i].nb_neighbors=k;
			}
		
		k=H->entry[j].nb_neighbors;
		ts=gonl_arra_govj(H->entry[j].neighbor,k,i,&dummy);
		if(ts==0)
			{if(k>=maxval)
				{fprintf(tmpout,"2-maxval=%d is reached\n",maxval);
				exit(0);
				}
			H->entry[j].neighbor[k]=i;
			k++;
			H->entry[j].nb_neighbors=k;
			}
		}
	}

for(i=0;i<N;i++)
	{w=0;
	for(j=0;j<N;j++)if(j!=i)
		{dis1=wodt_dist_gilq(S[i].zent,S[j].zent);
		dis2=S[i].rad+S[j].rad+2.0*probe;
		if(dis1<=dis2)
			{if(w>=max_inter)
				{fprintf(tmpout,"max_inter=%d is reached\n",max_inter);
				exit(0);
				}
			ts=gonl_arra_govj(H->entry[i].neighbor,H->entry[i].nb_neighbors,j,&dummy);
			if(ts==0)
				{H->inter[i].neighbor[w]=j;
				w++;
				}
			}
		}
	H->inter[i].nb_neighbors=w;
	}
}



