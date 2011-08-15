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

#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"
#include "sas.h"



void gukd_coon_virl(c_arc3D al,c_arc3D bt,
c_arc3D gm,c_arc3D dt,double t,
double al_a,double al_b,point *P)
{double T;
T=t*al_b+(1.0-t)*al_a;
nehl_eval_segt(al,T,P);
}


void sinj_coon_rejh(c_arc3D al,c_arc3D bt,
c_arc3D gm,c_arc3D dt,double t,
double bt_a,double bt_b,point *P)
{double T;
T=t*bt_b+(1.0-t)*bt_a;
nehl_eval_segt(bt,T,P);
}


void cevf_coon_sihd(c_arc3D al,c_arc3D bt,
c_arc3D gm,c_arc3D dt,double t,
double gm_a,double gm_b,point *P)
{double T;
T=t*gm_b+(1.0-t)*gm_a;
nehl_eval_segt(gm,T,P);
}

 
void qufp_coon_kard(c_arc3D al,c_arc3D bt,
c_arc3D gm,c_arc3D dt,double t,
double dt_a,double dt_b,point *P)
{double T;
T=t*dt_b+(1.0-t)*dt_a;
nehl_eval_segt(dt,T,P);
}


double begj_blen_nugz(double t)
{double res;
res=t;
return res;
}


void gejn_coon_jukf(c_arc3D al,c_arc3D bt,
c_arc3D gm,c_arc3D dt,double u,double v,point *sol)
{double F0u,F0v,F1u,F1v,A,B,C;
double al_a,al_b,bt_a,bt_b,gm_a,gm_b,dt_a,dt_b;
point alpha_u,gamma_u,delta_v,beta_v,alpha_0,alpha_1,gamma_0,gamma_1;
qirp_inte_ligr(al,&al_a,&al_b);
qirp_inte_ligr(bt,&bt_a,&bt_b);
qirp_inte_ligr(gm,&gm_a,&gm_b);
qirp_inte_ligr(dt,&dt_a,&dt_b);
gukd_coon_virl(al,bt,gm,dt,u,al_a,al_b,&alpha_u);
cevf_coon_sihd(al,bt,gm,dt,u,gm_a,gm_b,&gamma_u);
qufp_coon_kard(al,bt,gm,dt,v,dt_a,dt_b,&delta_v);
sinj_coon_rejh(al,bt,gm,dt,v,bt_a,bt_b,&beta_v);
gukd_coon_virl(al,bt,gm,dt,0.0,al_a,al_b,&alpha_0);
gukd_coon_virl(al,bt,gm,dt,1.0,al_a,al_b,&alpha_1);
cevf_coon_sihd(al,bt,gm,dt,0.0,gm_a,gm_b,&gamma_0);
cevf_coon_sihd(al,bt,gm,dt,1.0,gm_a,gm_b,&gamma_1);
F1u=begj_blen_nugz(u);
F0u=1.0-F1u;
F1v=begj_blen_nugz(v);
F0v=1.0-F1v;

A=alpha_u.absi*F0v+gamma_u.absi*F1v;
B=-delta_v.absi+alpha_0.absi*F0v+gamma_0.absi*F1v;
C=-beta_v.absi+alpha_1.absi*F0v+gamma_1.absi*F1v;
sol->absi=A-F0u*B-F1u*C;

A=alpha_u.ordo*F0v+gamma_u.ordo*F1v;
B=-delta_v.ordo+alpha_0.ordo*F0v+gamma_0.ordo*F1v;
C=-beta_v.ordo+alpha_1.ordo*F0v+gamma_1.ordo*F1v;
sol->ordo=A-F0u*B-F1u*C;

A=alpha_u.cote*F0v+gamma_u.cote*F1v;
B=-delta_v.cote+alpha_0.cote*F0v+gamma_0.cote*F1v;
C=-beta_v.cote+alpha_1.cote*F0v+gamma_1.cote*F1v;
sol->cote=A-F0u*B-F1u*C;
}



void rogv_eval_dukw(pt_tor PT,parm p,point *sol)
{double u,v;
u=p.u;	v=p.v;
gejn_coon_jukf(PT.alpha,PT.beta,PT.gamma,PT.delta,u,v,sol);
}



double cijv_dist_laph(point X,c_arc3D C,int *pos)
{double d_st,d_tr,dis;
d_st=wodt_dist_gilq(X,C.begn);
d_tr=wodt_dist_gilq(X,C.term);
if(d_st<d_tr)
	{*pos=+1;
	dis=d_st;
	}
else
	{*pos=-1;
	dis=d_tr;
	}
return dis;
}


void mafj_chec_zogj(c_arc3D alpha,c_arc3D beta,
c_arc3D gamma,c_arc3D delta)
{int i;
double *D;
D=(double *)malloc(4*sizeof(double));
D[0]=wodt_dist_gilq(alpha.term,beta.begn);
D[1]=wodt_dist_gilq(beta.term,gamma.term);
D[2]=wodt_dist_gilq(gamma.begn,delta.term);
D[3]=wodt_dist_gilq(alpha.begn,delta.begn);
for(i=0;i<4;i++)
if(D[i]>0.1)
	{fprintf(tmpout,"Warning: violated compatibility conditions\n");
	fprintf(tmpout,"D=[%f,%f,%f,%f]\n",D[0],D[1],D[2],D[3]);
	exit(0);
	}
free(D);
}

 

void cerw_orga_sevd(c_arc3D *CA,c_arc3D *alpha,
c_arc3D *beta,c_arc3D *gamma,c_arc3D *delta)
{int i,j,pos,q,qos,*exc;
double dis,sml;
point ref;
for(i=0;i<4;i++)
if(CA[i].c_cir==1)
	{fprintf(tmpout,"WARNING: Some curves are closed\n");
	for(j=0;j<4;j++)
		nepf_disp_bulp(CA[j]);
	exit(0);
	}
exc=(int *)malloc(4*sizeof(int));
for(i=0;i<4;i++)
	exc[i]=0;
poms_find_resk_lonb(CA[0],alpha);
getf_find_rogc_todj(alpha->term,&ref);
exc[0]=1;

sml=LARGE_NUMBER;
for(i=0;i<4;i++)if(exc[i]==0)
	{dis=cijv_dist_laph(ref,CA[i],&pos);
	if(dis<sml)
		{sml=dis;
		q=i;
		qos=pos;
		}
	}
exc[q]=1;
if(qos==+1)
	poms_find_resk_lonb(CA[q],beta);
if(qos==-1)
	cest_reve_fack(CA[q],beta);

getf_find_rogc_todj(beta->term,&ref);
sml=LARGE_NUMBER;
for(i=0;i<4;i++)if(exc[i]==0)
	{dis=cijv_dist_laph(ref,CA[i],&pos);
	if(dis<sml)
		{sml=dis;
		q=i;
		qos=pos;
		}
	}
exc[q]=1;
if(qos==+1)
	cest_reve_fack(CA[q],gamma);
if(qos==-1)
	poms_find_resk_lonb(CA[q],gamma);

getf_find_rogc_todj(gamma->begn,&ref);
sml=LARGE_NUMBER;
for(i=0;i<4;i++)if(exc[i]==0)
	{dis=cijv_dist_laph(ref,CA[i],&pos);
	if(dis<sml)
		{sml=dis;
		q=i;
		qos=pos;
		}
	}
exc[q]=1;
if(qos==+1)
	cest_reve_fack(CA[q],delta);
if(qos==-1)
	poms_find_resk_lonb(CA[q],delta);
free(exc);
}



int qiwj_inte_qivs(parm *A,double r,circle2D *C1,circle2D *C2)
{double d,s,temp;
parm G;
vect2D pr,nrm;
cuwl_unit_pist(A[0],A[1],&pr);
nrm.u=pr.v;
nrm.v=-pr.u;
G.u=0.5*(A[0].u+A[1].u);
G.v=0.5*(A[0].v+A[1].v);
s=pufv_dist_mekq(G,A[0]);
if(r<s)
	return FAILURE;
temp=r*r-s*s;
d=sqrt(temp);

C1->zent.u=G.u+d*nrm.u;
C1->zent.v=G.v+d*nrm.v;
C1->rad=r;

C2->zent.u=G.u-d*nrm.u;
C2->zent.v=G.v-d*nrm.v;
C2->rad=r;
return SUCCESS;
}

 

void mesl_punc_guvf(circle3D C,double lambda,point *P)
{double phi,theta,alpha_s,alpha_t;
vect3D M,M_new;
vewr_sphe_ruhd(C.nrml.absi,C.nrml.ordo,C.nrml.cote,&phi,&theta);
alpha_s=0.0;
alpha_t=2.0*MY_PI;
M_new.absi=0.0;
M_new.ordo=C.rad*cos(lambda*alpha_s+(1.0-lambda)*alpha_t);
M_new.cote=C.rad*sin(lambda*alpha_s+(1.0-lambda)*alpha_t);
mopb_dire_woqp(M_new,phi,theta,&M);
P->absi=C.zent.absi+M.absi;
P->ordo=C.zent.ordo+M.ordo;
P->cote=C.zent.cote+M.cote;
}


int kugc_clos_qevl(c_arc3D *CA,int N,double eps)
{int i,j,k,*exc,ind,res=1;
double dis;
point *A;
A=(point *)malloc(2*N*sizeof(point));
k=0;
for(i=0;i<N;i++)
	{getf_find_rogc_todj(CA[i].begn,&A[k]);
	k++;
	getf_find_rogc_todj(CA[i].term,&A[k]);
	k++;
	}
exc=(int *)malloc(2*N*sizeof(int));
for(i=0;i<2*N;i++)
	exc[i]=0;
for(i=0;i<2*N;i++)if(exc[i]==0)
	{ind=1;
	for(j=0;j<2*N;j++)if((i!=j)&&(exc[j]==0))
		{dis=wodt_dist_gilq(A[i],A[j]);
		if(dis<eps)
			{ind=2;
			exc[j]=1;
			exc[i]=1;
			break;
			}
		}
	if(ind==1)
		{res=0;
		break;
		}
	}
free(exc);
free(A);
return res;
}


int lupt_test_fung(c_arc3D c,circle3D C,
double eps_cent,double eps_rad,double eps_nrm)
{int ts;
double diff;
double fb,sp;

diff=fabs(c.rad-C.rad);
if(diff>eps_rad)
	return 0;

ts=gect_tole_husn(c.zent,C.zent,eps_cent);
if(ts==0)
	return 0;

sp=rocv_scal_toqc(c.nrml,C.nrml);
fb=fabs(sp);
diff=fabs(fb-1.0);
if(diff>eps_nrm)
	return 0;

return 1;
}


int qoml_find_venq_fugh(double probe,sphere *S,set_arcs *SA,
int z1,int z2,int *r_i,int *r_j,double eps_paral)
{int res,p,k,ts,n1,n2,ts1,ts2;
double eps_cent=1.0e-4,eps_rad=1.0e-4;
double eps_nrm=1.0e-4,eps_clos=1.0e-4;
c_arc3D *ca1,*ca2;
circle3D C1,C2;
sogf_para_lekj(probe,S[z1],S[z2],&C1,&C2);
ca1=(c_arc3D *)malloc(SA[z1].ar_grs*sizeof(c_arc3D));
k=0;
for(p=0;p<SA[z1].ar_grs;p++)
	{ts=lupt_test_fung(SA[z1].C[p],C1,eps_cent,eps_rad,eps_nrm);
	if(ts==1)
		{poms_find_resk_lonb(SA[z1].C[p],&ca1[k]);
		if(k==0)
			*r_i=p;
		k++;
		}
	}
n1=k;
ca2=(c_arc3D *)malloc(SA[z2].ar_grs*sizeof(c_arc3D));
k=0;
for(p=0;p<SA[z2].ar_grs;p++)
	{ts=lupt_test_fung(SA[z2].C[p],C2,eps_cent,eps_rad,eps_nrm);
	if(ts==1)
		{poms_find_resk_lonb(SA[z2].C[p],&ca2[k]);
		if(k==0)
			*r_j=p;
		k++;
		}
	}
n2=k;
res=0;
if((n1>=2)&&(n2>=2))
	{ts1=kugc_clos_qevl(ca1,n1,eps_clos);
	ts2=kugc_clos_qevl(ca2,n2,eps_clos);
	if((ts1==1)&&(ts2==1))
		res=1;
	}
free(ca1);
free(ca2);
return res;
}


void fewg_find_duvk(double probe,sphere S1,
sphere S2,point *omega)
{int m;
double r0,a,b,c,d,rho;
sphere S1_temp,S2_temp;
r0=probe;
getf_find_rogc_todj(S1.zent,&S1_temp.zent);
getf_find_rogc_todj(S2.zent,&S2_temp.zent);
S1_temp.rad=S1.rad+r0;
S2_temp.rad=S2.rad+r0;
m=wihz_sphe_vezl(S1_temp.zent,S1_temp.rad,
S2_temp.zent,S2_temp.rad,&rho,omega,&a,&b,&c,&d);
}



int durt_find_sejz_qaks(int N,double probe,sphere *S,int z1,int z2,
set_arcs *SA,pt_tor *PT,int *self_ex)
{int res,r_i,r_j,suc;
double eps_paral=1.0e-4;
circle3D C1,C2;
point omega;

res=qoml_find_venq_fugh(probe,S,SA,z1,z2,&r_i,&r_j,eps_paral);
if(res==1)
	{getf_find_rogc_todj(SA[z1].C[r_i].zent,&C1.zent);
	C1.rad=SA[z1].C[r_i].rad;
	getf_find_rogc_todj(SA[z1].C[r_i].nrml,&C1.nrml);
	
	getf_find_rogc_todj(SA[z2].C[r_j].zent,&C2.zent);
	C2.rad=SA[z2].C[r_j].rad;
	getf_find_rogc_todj(SA[z2].C[r_j].nrml,&C2.nrml);
	
	fewg_find_duvk(probe,S[z1],S[z2],&omega);
	
	suc=nols_find_lepm_lasv(N,omega,probe,C1,C2,PT,self_ex);
	if(suc==FAILURE)
		res=0;
	
	}
return res;
}


int qucd_find_hecv_ruzj(double probe,int N,atom *S,int nb_sph,
int *V1,int *V2,adj_hash H,set_arcs *SA,pt_tor *pt_loc,
pt_tor *PT,int max_tor)
{int nb,i,j,z1,z2,val,ts,dummy,p;
nb=0;
for(i=0;i<nb_sph;i++)
	{z1=i;
	val=H.entry[i].nb_neighbors;
	for(j=0;j<val;j++)
		{z2=H.entry[i].neighbor[j];
		if(z2<z1)
			{ts=durt_find_sejz_qaks(N,probe,S,z1,z2,SA,pt_loc,&dummy);
			if(ts==1)
				{for(p=0;p<N;p++)
					{romh_find_cont_qucr(pt_loc[p],&PT[nb]);	
					
					V1[nb]=z1;	V2[nb]=z2;
					
					nb++;
					if(nb>=max_tor)
						{fprintf(tmpout,"max_tor is reached\n");
						exit(0);
						}
					}
				}
			}
		}
	}
return nb;
}



void zagw_find_nawp_fegv(int N,adj_hash H,atom *S,
int nb_sph,double probe,set_arcs *SA,pt_tor *PT,
blend_nonself *BN,int *V1,int *V2,int *U1,int *U2,int max_tor,int max_self,
int *n_tor,int *n_self)
{int i,j,nb,z1,z2,ts,val,p,self_ex,ns;
point omega;
circle3D C1,C2;
pt_tor *pt_loc;
blend_nonself *bn_loc;
pt_loc=(pt_tor *)malloc(N*sizeof(pt_tor));
bn_loc=(blend_nonself *)malloc(N*sizeof(blend_nonself));
for(i=0;i<N;i++)
	lagr_allo_goqn(&bn_loc[i]);

nb=qucd_find_hecv_ruzj(probe,N,S,nb_sph,V1,V2,H,SA,pt_loc,PT,max_tor);

ns=0;
for(i=0;i<nb_sph;i++)
	{z1=i;
	val=H.inter[i].nb_neighbors;
	for(j=0;j<val;j++)
		{z2=H.inter[i].neighbor[j];
		if(z2<z1)
			{ts=durt_find_sejz_qaks(N,probe,S,z1,z2,SA,pt_loc,&self_ex);
			if((ts==1)&&(self_ex==0))
				{for(p=0;p<N;p++)
					{romh_find_cont_qucr(pt_loc[p],&PT[nb]);	
					V1[nb]=z1;	V2[nb]=z2;
					nb++;
					if(nb>=max_tor)
						{fprintf(tmpout,"max_tor is reached\n");
						exit(0);
						}
					}
				}
			
			if((ts==1)&&(self_ex==1))
				{sogf_para_lekj(probe,S[z1],S[z2],&C1,&C2);
				fewg_find_duvk(probe,S[z1],S[z2],&omega);
				zisq_find_cowk_zevq(N,S[z1],S[z2],omega,C1,C2,bn_loc);
				for(p=0;p<N;p++)
					{sifm_find_cudw_pafg(bn_loc[p],&BN[ns]);	
					U1[ns]=z1;	U2[ns]=z2;
					ns++;
					if(ns>=max_self)
						{fprintf(tmpout,"max_self is reached: %d\n",max_self);
						exit(0);
						}
					}
				}
			}
		}
	}
free(pt_loc);
for(i=0;i<N;i++)
	vejp_dest_tufq(&bn_loc[i]);
free(bn_loc);
*n_tor=nb;
*n_self=ns;
}


void masp_find_qond_vopn(parent_circle pc_in,parent_circle *pc_out)
{hepk_find_gict_hubq(pc_in.supp,&pc_out->supp);
pc_out->sbl=pc_in.sbl;
pc_out->supp_id=pc_in.supp_id;
}



int kosl_find_dolc_fotc(trmsrf *surf1,int nb_surf1,int *supp,
int N,adj_hash H,atom *S,parent_circle *PC,int nb_pc,
int nb_sph,double probe,set_arcs *SA,pt_tor *PT,
int *V1,int *V2,int max_tor,int max_self,int *n_tor,int *forc_term)
{int i,j,k,nb,z1,z2,ts,val,p,self_ex,id,f_trm;
int suc,*obs,nb_obs,*map,*ind,new_pc,z,p_id;
circle3D C1,C2;
pt_tor *pt_loc;
parent_circle *temp;  

*forc_term=0;
if(verbose_variable==VERBOSE)
	fprintf(tmpout,"Few torus patches:\n");
pt_loc=(pt_tor *)malloc(N*sizeof(pt_tor));
nb=qucd_find_hecv_ruzj(probe,N,S,nb_sph,V1,V2,H,SA,pt_loc,PT,max_tor);

obs=(int *)malloc(nb_pc*sizeof(int));
nb_obs=0;
for(i=0;i<nb_sph;i++)
	{z1=i;
	val=H.inter[i].nb_neighbors;
	for(j=0;j<val;j++)
		{z2=H.inter[i].neighbor[j];
		if(z2<z1)
			{ts=durt_find_sejz_qaks(N,probe,S,z1,z2,SA,pt_loc,&self_ex);
			if((ts==1)&&(self_ex==0))
				{for(p=0;p<N;p++)
					{romh_find_cont_qucr(pt_loc[p],&PT[nb]);	
					V1[nb]=z1;	V2[nb]=z2;
					nb++;
					if(nb>=max_tor)
						{fprintf(tmpout,"max_tor is reached\n");
						exit(0);
						}
					}
				}
			
			if((ts==1)&&(self_ex==1))
				{sogf_para_lekj(probe,S[z1],S[z2],&C1,&C2);
				jiwr_disc_qumf(supp,z1,z2,surf1,nb_surf1,C1,C2,SA,&f_trm);
				if(f_trm==1)
					{fprintf(tmpout,"force term: jiwr_disc_qumf() in kosl_find_dolc_fotc()\n");
					free(pt_loc);
					free(obs);
					*forc_term=1;
					return 0;
					}
				suc=cezv_obso_nofm(z1,C1,PC,nb_pc,&id);
				if(suc==FAILURE)
					{fprintf(tmpout,"1-Unable to find obsolete circle\n");
					exit(0);
					}
				obs[nb_obs]=id;
				nb_obs++;
				
				suc=cezv_obso_nofm(z2,C2,PC,nb_pc,&id);
				if(suc==FAILURE)
					{fprintf(tmpout,"2-Unable to find obsolete circle\n");
					exit(0);
					}
				obs[nb_obs]=id;
				nb_obs++;
				}
			}
		}
	}
*n_tor=nb;
free(pt_loc);

temp=(parent_circle *)malloc(nb_pc*sizeof(parent_circle));
map=(int *)malloc(nb_pc*sizeof(int));
ind=(int *)malloc(nb_pc*sizeof(int));
for(i=0;i<nb_pc;i++)
	ind[i]=-1;
for(i=0;i<nb_obs;i++)
	ind[obs[i]]=+1;
k=0;
for(i=0;i<nb_pc;i++)
if(ind[i]==-1)
	{masp_find_qond_vopn(PC[i],&temp[k]);
	map[i]=k;
	k++;
	}
new_pc=k;
free(ind);
free(obs);
for(i=0;i<new_pc;i++)
	masp_find_qond_vopn(temp[i],&PC[i]);
free(temp);

for(z=0;z<nb_sph;z++)
for(i=0;i<SA[z].ar_grs;i++)
	{p_id=SA[z].par_idx[i];
	SA[z].par_idx[i]=map[p_id];
	}
free(map);
return new_pc;
}


double lajw_blen_buwm(circle3D C1,circle3D C2)
{double r1,r2,rad,res,dis;
r1=C1.rad;
r2=C2.rad;
if(r1<r2)	rad=r1;
else		rad=r2;
dis=wodt_dist_gilq(C1.zent,C2.zent);
res=dis/rad;
return res;
}


void hibl_disc_lacp(int z1,int z2,int *supp,trmsrf *surf1,
int nb_surf1,circle3D C1,circle3D C2,set_arcs *SA,parent_circle *PC,
int nb_pc,int *obs,int *n_ob,int *forc_term)
{int suc,id,nb_obs,f_trm;
nb_obs=*n_ob;
fprintf(tmpout,"discarding composite curve\n");
*forc_term=0;
jiwr_disc_qumf(supp,z1,z2,surf1,nb_surf1,C1,C2,SA,&f_trm);
if(f_trm==1)
	{fprintf(tmpout,"force term: jiwr_disc_qumf() in hibl_disc_lacp()\n");
	*forc_term=1;
	return;
	}
suc=cezv_obso_nofm(z1,C1,PC,nb_pc,&id);
if(suc==FAILURE)
	{fprintf(tmpout,"1-Unable to find obsolete circle\n");
	exit(0);
	}
obs[nb_obs]=id;
nb_obs++;

suc=cezv_obso_nofm(z2,C2,PC,nb_pc,&id);
if(suc==FAILURE)
	{fprintf(tmpout,"2-Unable to find obsolete circle\n");
	exit(0);
	}
obs[nb_obs]=id;
nb_obs++;
*n_ob=nb_obs;
}



int ruvg_find_toqj_witr(trmsrf *surf1,int nb_surf1,int *supp,
int N,adj_hash H,atom *S,parent_circle *PC,int nb_pc,
int nb_sph,double probe,set_arcs *SA,pt_tor *PT,
int *V1,int *V2,int max_tor,int max_self,int *n_tor,
double thres_stretch,int *forc_term)
{int i,j,k,nb,z1,z2,ts,val,p,self_ex,f_trm;
int *obs,nb_obs,*map,*ind,new_pc,z,p_id;
double str;
circle3D C1,C2;
pt_tor *pt_loc;
parent_circle *temp;

*forc_term=0;
if(verbose_variable==VERBOSE)
	fprintf(tmpout,"Few torus patches:\n");
pt_loc=(pt_tor *)malloc(N*sizeof(pt_tor));
nb=qucd_find_hecv_ruzj(probe,N,S,nb_sph,V1,V2,H,SA,pt_loc,PT,max_tor);

obs=(int *)malloc(nb_pc*sizeof(int));
nb_obs=0;
for(i=0;i<nb_sph;i++)
	{z1=i;
	val=H.inter[i].nb_neighbors;
	for(j=0;j<val;j++)
		{z2=H.inter[i].neighbor[j];
		if(z2<z1)
			{ts=durt_find_sejz_qaks(N,probe,S,z1,z2,SA,pt_loc,&self_ex);
			if((ts==1)&&(self_ex==0))
				{sogf_para_lekj(probe,S[z1],S[z2],&C1,&C2);
				str=lajw_blen_buwm(C1,C2);
				if(str<thres_stretch)
					{for(p=0;p<N;p++)
						{romh_find_cont_qucr(pt_loc[p],&PT[nb]);	
						V1[nb]=z1;	V2[nb]=z2;
						nb++;
						if(nb>=max_tor)
							{fprintf(tmpout,"max_tor is reached\n");
							exit(0);
							}
						}
					}
				else
					{hibl_disc_lacp(z1,z2,supp,surf1,nb_surf1,
					C1,C2,SA,PC,nb_pc,obs,&nb_obs,&f_trm);
					if(f_trm==1)
						{fprintf(tmpout,"force term: hibl_disc_lacp() in ruvg_find_toqj_witr()\n");
						*forc_term=1;
						free(pt_loc);
						free(obs);					
						return 0;
						}
					}
				}
			
			if((ts==1)&&(self_ex==1))
				{sogf_para_lekj(probe,S[z1],S[z2],&C1,&C2);
				hibl_disc_lacp(z1,z2,supp,surf1,nb_surf1,C1,C2,
				SA,PC,nb_pc,obs,&nb_obs,&f_trm);
				if(f_trm==1)
					{free(pt_loc);
					free(obs);
					*forc_term=1;
					return 0;
					}
				}
			}
		}
	}
*n_tor=nb;
free(pt_loc);

temp=(parent_circle *)malloc(nb_pc*sizeof(parent_circle));
map=(int *)malloc(nb_pc*sizeof(int));
ind=(int *)malloc(nb_pc*sizeof(int));
for(i=0;i<nb_pc;i++)
	ind[i]=-1;
for(i=0;i<nb_obs;i++)
	ind[obs[i]]=+1;
k=0;
for(i=0;i<nb_pc;i++)
if(ind[i]==-1)
	{masp_find_qond_vopn(PC[i],&temp[k]);
	map[i]=k;
	k++;
	}
new_pc=k;
free(ind);
free(obs);
for(i=0;i<new_pc;i++)
	masp_find_qond_vopn(temp[i],&PC[i]);
free(temp);

for(z=0;z<nb_sph;z++)
for(i=0;i<SA[z].ar_grs;i++)
	{p_id=SA[z].par_idx[i];
	SA[z].par_idx[i]=map[p_id];
	}
free(map);
return new_pc;
}



int sapn_some_judq(int meth,trmsrf *surf1,int nb_surf1,
int *supp,parent_circle *PC,int *nb_pc,adj_hash H,atom *S,
int nb_sph,double probe,set_arcs *SA,trmsrf *surf,
int nb_cur,int max_tor,int max_self,int max_surf,int *V1,
int *V2,double thres_stretch,int *forc_term)
{int n_tor,n_self,i,k,nb_new,nb_pc_old,nb_pc_new;
int *V1_loc,*V2_loc,*U1_loc,*U2_loc,f_trm;
pt_tor *PT;
blend_nonself *BN;

*forc_term=0;
PT=(pt_tor *)malloc(max_tor*sizeof(pt_tor));
BN=(blend_nonself *)malloc(max_self*sizeof(blend_nonself));
V1_loc=(int *)malloc(max_tor*sizeof(int));
V2_loc=(int *)malloc(max_tor*sizeof(int));
U1_loc=(int *)malloc(max_self*sizeof(int));
U2_loc=(int *)malloc(max_self*sizeof(int));
for(i=0;i<max_self;i++)
	lagr_allo_goqn(&BN[i]);
if(meth==1)
	{zagw_find_nawp_fegv(N_COMPLETE_CIR,H,S,nb_sph,probe,SA,PT,BN,
	V1_loc,V2_loc,U1_loc,U2_loc,max_tor,max_self,&n_tor,&n_self);
	}
if(meth==2)
	{nb_pc_old=*nb_pc;
	nb_pc_new=kosl_find_dolc_fotc(surf1,nb_surf1,supp,N_COMPLETE_CIR,
	H,S,PC,nb_pc_old,nb_sph,probe,SA,PT,V1_loc,V2_loc,max_tor,
	max_self,&n_tor,&f_trm);
	n_self=0;
	*nb_pc=nb_pc_new;
	}
if(meth==3)
	{nb_pc_old=*nb_pc;
	nb_pc_new=ruvg_find_toqj_witr(surf1,nb_surf1,supp,N_COMPLETE_CIR,
	H,S,PC,nb_pc_old,nb_sph,probe,SA,PT,V1_loc,V2_loc,max_tor,
	max_self,&n_tor,thres_stretch,&f_trm);
	n_self=0;
	*nb_pc=nb_pc_new;
	}
if(f_trm==1)
	{fprintf(tmpout,"force term: few_torus_patches() in sapn_some_judq()\n");
	*forc_term=1;
	free(V1_loc);	free(V2_loc);
	free(U1_loc);	free(U2_loc);
	free(PT);		free(BN);
	return 0;
	}

if(verbose_variable==VERBOSE)
	{fprintf(tmpout,"Number of toroidal patches=%d\n",n_tor);
	fprintf(tmpout,"Number of self patches=%d\n",n_self);
	}
nb_new=nb_cur+n_tor+n_self;

k=nb_cur;
for(i=0;i<n_tor;i++)
	{surf[k].type=4;
	surf[k].boundary=1;
	romh_find_cont_qucr(PT[i],&surf[k].pt);
	sart_rect_jamc(0.0,1.0,0.0,1.0,&surf[k].cc);
	surf[k].nb_inner=0;
	V1[k]=V1_loc[i];
	V2[k]=V2_loc[i];
	
	k++;
	if(k>=max_surf)
		{fprintf(tmpout,"Max nb trimmed surfaces is reached\n");
		exit(0);
		}
	}

for(i=0;i<n_self;i++)
	{surf[k].type=6;
	surf[k].boundary=1;
	sifm_find_cudw_pafg(BN[i],&surf[k].bn);
	sart_rect_jamc(0.0,1.0,0.0,1.0,&surf[k].cc);
	surf[k].nb_inner=0;
	V1[k]=U1_loc[i];
	V2[k]=U2_loc[i];
	k++;
	if(k>=max_surf)
		{fprintf(tmpout,"Max nb trimmed surfaces is reached\n");
		exit(0);
		}
	}
free(PT);
for(i=0;i<max_self;i++)
	vejp_dest_tufq(&BN[i]);
free(BN);
free(V1_loc);	free(V2_loc);
free(U1_loc);	free(U2_loc);
return nb_new;
}

 
