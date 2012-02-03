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


double vuqg_dist_faql(pt_tor PT,sphere S)
{int i;
double dis,res;
c_arc3D *CA;
CA=(c_arc3D *)malloc(4*sizeof(c_arc3D));
poms_find_resk_lonb(PT.alpha,&CA[0]);
poms_find_resk_lonb(PT.beta,&CA[1]);
poms_find_resk_lonb(PT.gamma,&CA[2]);
poms_find_resk_lonb(PT.delta,&CA[3]);
res=LARGE_NUMBER;
for(i=0;i<4;i++)
	{dis=zadt_erro_vehp(S,CA[i]);
	if(dis<res)
		res=dis;
	}
free(CA);
return res;
}


void vupn_chec_wudq(trmsrf *surf,
sphere *S,blend_tor BT)
{int r[2],i,z;
double dis,eps=1.0e-3;
r[0]=BT.sph_idx1;
r[1]=BT.sph_idx2;
z=BT.trim_idx;
for(i=0;i<2;i++)
	{dis=vuqg_dist_faql(surf[z].pt,S[r[i]]);
	if(dis>eps)
		{fprintf(tmpout,"dis=%f\n",dis);
		fprintf(tmpout,"side atoms [%d,%d]\n",r[0],r[1]);
		fprintf(tmpout,"Too large blend interface\n");
		exit(0);
		}
	}
}



void sedr_fill_luch(int nb_sph,blend_cpx *BC)
{int i,j,w,ts,dummy,v1,v2;
for(i=0;i<nb_sph;i++)
	BC->HE[i].nb=0;
for(j=0;j<BC->bt_grs;j++)
	{v1=BC->BT[j].sph_idx1;
	w=BC->HE[v1].nb;
	if(w>=MAX_LOC_BLEND)
		{fprintf(tmpout,"Current number of blend tor=%d\n",BC->bt_grs);
		fprintf(tmpout,"MAX_LOC_BLEND is exceeded\n");
		exit(0);
		}
	ts=gonl_arra_govj(BC->HE[v1].list,BC->HE[v1].nb,j,&dummy);
	if(ts==0)
		{BC->HE[v1].list[w]=j;
		BC->HE[v1].nb=w+1;
		}
	
	v2=BC->BT[j].sph_idx2;
	w=BC->HE[v2].nb;
	if(w>=MAX_LOC_BLEND)
		{fprintf(tmpout,"MAX_LOC_BLEND is exceeded\n");
		exit(0);
		}
	ts=gonl_arra_govj(BC->HE[v2].list,BC->HE[v2].nb,j,&dummy);
	if(ts==0)
		{BC->HE[v2].list[w]=j;
		BC->HE[v2].nb=w+1;
		}
	}
fprintf(tmpout,"Filling is complete\n");
}



void kejh_fill_zogq(sphere *S,trmsrf *surf2,
int nb_sph,blend_cpx *BC)
{int i,z;
fprintf(tmpout,"Filling blend complex\n");
for(i=0;i<BC->bt_grs;i++)
	{vupn_chec_wudq(surf2,S,BC->BT[i]);
	z=BC->BT[i].trim_idx;
	mafj_chec_zogj(surf2[z].pt.alpha,surf2[z].pt.beta,
	surf2[z].pt.gamma,surf2[z].pt.delta);
	}
if(BC->bt_grs>=MAX_BLEND_TOR)
	{fprintf(tmpout,"MAX_BLEND_TOR has been exceeded\n");
	exit(0);
	}
fprintf(tmpout,"Good blend integrity\n");
sedr_fill_luch(nb_sph,BC);
}



void rows_eval_qusg(pt_cnv PC,parm p,point *sol)
{double mu;
point val1,val2,A,B,C;
baryc2D q;
q.lambda1=1.0-p.u-p.v;
q.lambda2=p.u;
q.lambda3=p.v;
if(PC.sitn==0)
	suzr_find_heqz_wikr(q,PC.alpha,PC.beta,PC.gamma,sol);
if(PC.sitn==1)
	{suzr_find_heqz_wikr(q,PC.alpha,PC.beta,PC.gamma,&val1);
	getf_find_rogc_todj(PC.alpha.begn,&A);
	getf_find_rogc_todj(PC.beta.begn ,&B);
	getf_find_rogc_todj(PC.gamma.begn,&C);
	val2.absi=q.lambda1*A.absi+q.lambda2*B.absi+q.lambda3*C.absi;
	val2.ordo=q.lambda1*A.ordo+q.lambda2*B.ordo+q.lambda3*C.ordo;
	val2.cote=q.lambda1*A.cote+q.lambda2*B.cote+q.lambda3*C.cote;
	mu=q.lambda1*q.lambda2*q.lambda3;
	mu=PC.scl*mu;
	sol->absi=(1.0-mu)*val1.absi+mu*val2.absi;
	sol->ordo=(1.0-mu)*val1.ordo+mu*val2.ordo;
	sol->cote=(1.0-mu)*val1.cote+mu*val2.cote;
	}
if(PC.sitn==2)
	dufj_eval_wejf(PC.T,p.u,p.v,sol);
}


void vesj_reve_jekc(c_arc3D *alpha,c_arc3D *beta,
c_arc3D *gamma)
{c_arc3D al,bt,gm;
cest_reve_fack(*alpha,&al);
cest_reve_fack(*gamma,&bt);
cest_reve_fack(*beta ,&gm);
poms_find_resk_lonb(al,alpha);
poms_find_resk_lonb(bt,beta);
poms_find_resk_lonb(gm,gamma);
}


void qewd_rege_leht(sphere S,point *X)
{vect3D U;
culm_unit_peks(S.zent,*X,&U);
X->absi=S.zent.absi+S.rad*U.absi;
X->ordo=S.zent.ordo+S.rad*U.ordo;
X->cote=S.zent.cote+S.rad*U.cote;
}


void rokl_refi_leqf(sphere S,c_arc3D al,c_arc3D bt,c_arc3D gm,
c_arc3D *al_son,c_arc3D *bt_son,c_arc3D *gm_son)
{point al_m,bt_m,gm_m;
renw_midp_mocw(al,&al_m);
renw_midp_mocw(bt,&bt_m);
renw_midp_mocw(gm,&gm_m);
qewd_rege_leht(S,&al_m);
qewd_rege_leht(S,&bt_m);
qewd_rege_leht(S,&gm_m);
dels_geod_tuzd(S.zent,S.rad,al_m,bt_m,al_son);
dels_geod_tuzd(S.zent,S.rad,bt_m,gm_m,bt_son);
dels_geod_tuzd(S.zent,S.rad,gm_m,al_m,gm_son);
}



void qovm_find_cihn_rifq(sphere S,int N,c_arc3D *alpha,
c_arc3D *beta,c_arc3D *gamma)
{int i;
c_arc3D al_par,al_son,bt_par,bt_son,gm_par,gm_son;
poms_find_resk_lonb(*alpha,&al_par);
poms_find_resk_lonb(*beta ,&bt_par);
poms_find_resk_lonb(*gamma,&gm_par);
for(i=0;i<N;i++)
	{rokl_refi_leqf(S,al_par,bt_par,gm_par,&al_son,&bt_son,&gm_son);
	poms_find_resk_lonb(al_son,&al_par);
	poms_find_resk_lonb(bt_son,&bt_par);
	poms_find_resk_lonb(gm_son,&gm_par);
	}
poms_find_resk_lonb(al_par,alpha);
poms_find_resk_lonb(bt_par,beta);
poms_find_resk_lonb(gm_par,gamma);
}



int melq_orie_gekw(double probe,c_arc3D alpha,
c_arc3D beta,c_arc3D gamma)
{int i,res,nb_shrink=5;
double err1,err2,sp;
point *A,*MD,G;
sphere S1,S2,S;
c_arc3D *CA;
vect3D W1,W2;

CA=(c_arc3D *)malloc(3*sizeof(c_arc3D));
A=(point *)malloc(3*sizeof(point));
MD=(point *)malloc(3*sizeof(point));
poms_find_resk_lonb(alpha,&CA[0]);
poms_find_resk_lonb(beta ,&CA[1]);
poms_find_resk_lonb(gamma,&CA[2]);
for(i=0;i<3;i++)
	{getf_find_rogc_todj(CA[i].begn,&A[i]);	
	renw_midp_mocw(CA[i],&MD[i]);
	}
nefr_inte_sujc(A,probe,&S1,&S2);
err1=0.0;	err2=0.0;
for(i=0;i<3;i++)
	{err1=err1+mulh_erro_cedm(S1,MD[i]);
	err2=err2+mulh_erro_cedm(S2,MD[i]);
	}
if(err1<err2)  neqg_find_lodr_bogm(S1,&S);
else		   neqg_find_lodr_bogm(S2,&S);

qovm_find_cihn_rifq(S,nb_shrink,&CA[0],&CA[1],&CA[2]);
for(i=0;i<3;i++)
	{getf_find_rogc_todj(CA[i].begn,&A[i]);	
	renw_midp_mocw(CA[i],&MD[i]);
	}
gotq_norm_bitg(A[0],A[1],A[2],&W1);
G.absi=(MD[0].absi+MD[1].absi+MD[2].absi)/3.0;
G.ordo=(MD[0].ordo+MD[1].ordo+MD[2].ordo)/3.0;
G.cote=(MD[0].cote+MD[1].cote+MD[2].cote)/3.0;
bofp_form_nukv(S.zent,G,&W2);
sp=rocv_scal_toqc(W1,W2);
res=+1;
if(sp<0.0)
	res=-1;
free(MD);
free(CA);
free(A);
return res;
}


void suwf_orga_cezm(double probe,c_arc3D *CA,
c_arc3D *alpha,c_arc3D *beta,c_arc3D *gamma)
{int i,pos,q,qos,*exc,ort;
double dis,sml;
point ref;
for(i=0;i<3;i++)
if(CA[i].c_cir==1)
	{fprintf(tmpout,"WARNING: Some curves are closed\n");
	exit(0);
	}
exc=(int *)malloc(3*sizeof(int));
for(i=0;i<3;i++)
	exc[i]=0;
poms_find_resk_lonb(CA[0],alpha);
getf_find_rogc_todj(alpha->term,&ref);
exc[0]=1;

sml=LARGE_NUMBER;
for(i=0;i<3;i++)if(exc[i]==0)
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

for(i=0;i<3;i++)if(exc[i]==0)
	{dis=cijv_dist_laph(ref,CA[i],&pos);
	sml=dis;
	q=i;
	qos=pos;
	break;
	}

if(qos==+1)
	poms_find_resk_lonb(CA[q],gamma);
if(qos==-1)
	cest_reve_fack(CA[q],gamma);
free(exc);

ort=melq_orie_gekw(probe,*alpha,*beta,*gamma);
if(ort==-1)
	vesj_reve_jekc(alpha,beta,gamma);
}


void jofd_find_mikn_gehj(pt_cnv P,pt_cnv *Q)
{poms_find_resk_lonb(P.alpha,&Q->alpha);
poms_find_resk_lonb(P.beta,&Q->beta);
poms_find_resk_lonb(P.gamma,&Q->gamma);
Q->sitn=P.sitn;
if(P.sitn==1)
	Q->scl=P.scl;
if(P.sitn==2)
	{zikt_find_jotz_jewb(P.T,&Q->T);
	kotg_find_wuhk_kemt(P.cc,&Q->cc);
	}
}


void qerw_corn_devf(trmsrf surf,point *A)
{getf_find_rogc_todj(surf.pc.alpha.begn,&A[0]);
getf_find_rogc_todj(surf.pc.beta.begn,&A[1]);
getf_find_rogc_todj(surf.pc.gamma.begn,&A[2]);
}


int cevw_amon_qutc(trmsrf *surf,int nb,double eps,point *corner)
{int res,i,j,dummy,ind,ts;
point *A;
A=(point *)malloc(3*sizeof(point));
res=0;
for(i=0;i<nb;i++)
	{qerw_corn_devf(surf[i],A);
	ind=1;
	for(j=0;j<3;j++)
		{ts=qidk_arra_ticg(A,3,corner[j],eps,&dummy);
		if(ts==0)
			{ind=2;
			break;
			}
		}
	if(ind==1)
		{res=1;
		break;
		}
	}
free(A);
return res;
}



int hicj_fill_zalj(double probe,sphere *S,int nb_sph,
trmsrf *surf2,int nb_cur,blend_cpx BC,
trmsrf *surf3,supp_sph_tri *treb,int max_surf)
{int v1,nb_loc,val,i,k,ts,n_allc;
double eps=1.0e-4;
point *corner;
c_arc3D *alpha,*beta,*gamma,*CA;
pt_cnv temp;
supp_sph_tri *treb_loc;

CA=(c_arc3D *)malloc(3*sizeof(c_arc3D));
corner=(point *)malloc(3*sizeof(point));
k=nb_cur;
for(v1=0;v1<nb_sph;v1++)
	{val=BC.HE[v1].nb;
	if(val==0)	n_allc=16;
	else		n_allc=val*val*16;
	alpha   = (c_arc3D *)malloc(n_allc*sizeof(c_arc3D));	
	beta    = (c_arc3D *)malloc(n_allc*sizeof(c_arc3D));	
	gamma   = (c_arc3D *)malloc(n_allc*sizeof(c_arc3D));	
	treb_loc= (supp_sph_tri *)malloc(n_allc*sizeof(supp_sph_tri));
	nb_loc=gazq_loca_joth(S,surf2,v1,BC,alpha,beta,gamma,treb_loc);
	if(verbose_variable==VERBOSE)
		{if(nb_loc!=0)
			fprintf(tmpout,"concave probe atom:   v1=%d   nb_loc=%d  nb_trim=%d\n",v1,nb_loc,k);
		}
	for(i=0;i<nb_loc;i++)
		{poms_find_resk_lonb(alpha[i],&CA[0]);
		poms_find_resk_lonb(beta[i]  ,&CA[1]);
		poms_find_resk_lonb(gamma[i] ,&CA[2]);
		suwf_orga_cezm(probe,CA,&temp.alpha,&temp.beta,&temp.gamma);
		getf_find_rogc_todj(temp.alpha.begn,&corner[0]);
		getf_find_rogc_todj(temp.beta.begn ,&corner[1]);
		getf_find_rogc_todj(temp.gamma.begn,&corner[2]);
		
		ts=cevw_amon_qutc(surf3,k,eps,corner);
		if(ts==0)
			{tunp_unit_noqr(&surf3[k].cc);
			jofd_find_mikn_gehj(temp,&surf3[k].pc);
			surf3[k].type=5;
			surf3[k].boundary=1;
			surf3[k].nb_inner=0;
			surf3[k].pc.sitn =0;
			treb[k].sph_idx1=treb_loc[i].sph_idx1;
			treb[k].sph_idx2=treb_loc[i].sph_idx2;
			treb[k].sph_idx3=treb_loc[i].sph_idx3;
			k++;
			if(k>=max_surf)
				{fprintf(tmpout,"Max nb trimmed surfaces is reached\n");
				exit(0);
				}
			}
		}
	free(treb_loc);
	free(alpha);
	free(beta);
	free(gamma);
	}
free(corner);
free(CA);
return k;
}


