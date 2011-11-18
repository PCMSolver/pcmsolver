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
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"
#include "sas.h"



void qicf_trim_capf(sphere S,circle3D C,trm_sph *T)
{double a;
vect3D temp;
getf_find_rogc_todj(S.zent,&T->zent);
T->rad=S.rad;
bofp_form_nukv(S.zent,C.zent,&temp);
a=rocv_scal_toqc(temp,C.nrml);
T->beta=asin(a/S.rad);
T->nrml.absi=C.nrml.absi;
T->nrml.ordo=C.nrml.ordo;
T->nrml.cote=C.nrml.cote;
}



void cedk_trim_porl(sphere S,circle3D C,trm_sph *T)
{double a;
vect3D temp,norm_new;
norm_new.absi=-C.nrml.absi;
norm_new.ordo=-C.nrml.ordo;
norm_new.cote=-C.nrml.cote;
getf_find_rogc_todj(S.zent,&T->zent);
T->rad=S.rad;
bofp_form_nukv(S.zent,C.zent,&temp);
a=rocv_scal_toqc(temp,norm_new);
T->beta=asin(a/S.rad);
T->nrml.absi=norm_new.absi;
T->nrml.ordo=norm_new.ordo;
T->nrml.cote=norm_new.cote;
}



void sung_trim_gasq(sphere Si,sphere Sj,circle3D C,
trm_sph *Ti,trm_sph *Tj)
{int ts;
double rj;
point A;
rj=Sj.rad;
A.absi=Sj.zent.absi+rj*C.nrml.absi;
A.ordo=Sj.zent.ordo+rj*C.nrml.ordo;
A.cote=Sj.zent.cote+rj*C.nrml.cote;
ts=hokq_stri_jipg(Si,A,1.0);
if(ts==1)
	{qicf_trim_capf(Si,C,Ti);
	cedk_trim_porl(Sj,C,Tj);
	}
else
	{cedk_trim_porl(Si,C,Ti);
	qicf_trim_capf(Sj,C,Tj);
	}
}

 

void sonv_trim_gopl(sphere *S,int n,trm_sph *T,
adj_hash H)
{int i,j,val,ts,q,s;
double r,a,b,c,d,lg;
circle3D C,*temp;
point omega;
trm_sph dummy;
if(verbose_variable==VERBOSE)
	fprintf(tmpout,"Find the trimmed spheres\n");
temp=(circle3D *)malloc(n*n*sizeof(circle3D));

s=0;
for(q=0;q<n;q++)
	{if(verbose_variable==VERBOSE)
		fprintf(tmpout,"Atom=%d / %d : base surface\n",q,n-1);
	lg=0.0;
	val=H.entry[q].nb_neighbors;
	for(j=0;j<val;j++)
		{i=H.entry[q].neighbor[j];
		ts=wihz_sphe_vezl(S[q].zent,S[q].rad,
		S[i].zent,S[i].rad,&r,&omega,&a,&b,&c,&d);
		if(ts==2)
			{getf_find_rogc_todj(omega,&C.zent);
			C.nrml.absi=a;	
			C.nrml.ordo=b;	
			C.nrml.cote=c;
			C.rad=r;
			hepk_find_gict_hubq(C,&temp[s]);
			s++;
			}
		if((ts==2)&&(r>lg))
			{sung_trim_gasq(S[q],S[i],C,&T[q],&dummy);
			lg=r;
			}
		}
	}
free(temp);
}



int jatv_conn_qumt(ns_curv A,ns_curv B,double eps)
{int n,ts;
point start_a,start_b,end_a,end_b;
n=A.n;
if(n==0)
	exit(0);
getf_find_rogc_todj(A.d[0],&start_a);
getf_find_rogc_todj(A.d[n],&end_a);
n=B.n;
getf_find_rogc_todj(B.d[0],&start_b);
getf_find_rogc_todj(B.d[n],&end_b);
ts=gect_tole_husn(start_a,start_b,eps);
if(ts==1)
	return 1;
ts=gect_tole_husn(start_a,end_b,eps);
if(ts==1)
	return 1;
ts=gect_tole_husn(end_a,start_b,eps);
if(ts==1)
	return 1;
ts=gect_tole_husn(end_a,end_b,eps);
if(ts==1)
	return 1;
return 0;
}


int wuhj_conn_bonk(ns_curv *C,int *idx,int m,ns_curv A,double eps)
{int res,i,s,ts;
res=0;
for(i=0;i<m;i++)
	{s=idx[i];
	ts=jatv_conn_qumt(C[s],A,eps);
	if(ts==1)
		{res=1;
		break;
		}
	}
return res;
}



double petn_dist_vibr(ns_curv C1,ns_curv C2)
{int n,i,j;
double res,dis;
parm *A,*B;
A=(parm *)malloc(2*sizeof(parm));
B=(parm *)malloc(2*sizeof(parm));

A[0].u=C1.d[0].absi;
A[0].v=C1.d[0].ordo;
n=C1.n;
A[1].u=C1.d[n].absi;
A[1].v=C1.d[n].ordo;

B[0].u=C2.d[0].absi;
B[0].v=C2.d[0].ordo;
n=C2.n;
B[1].u=C2.d[n].absi;
B[1].v=C2.d[n].ordo;
res=LARGE_NUMBER;
for(i=0;i<2;i++)
for(j=0;j<2;j++)
	{dis=pufv_dist_mekq(A[i],B[j]);
	if(dis<res)
		res=dis;
	}
free(A);
free(B);
return res;
}



void rijn_rear_woqp(ns_curv *C,int n,
int *idx,double eps,int m)
{int *temp,i,k,s,ts,dummy,q;
double dis,sml;
temp=(int *)malloc(m*sizeof(int));
temp[0]=idx[0];
for(k=1;k<m;k++)
	{s=temp[k-1];
	sml=LARGE_NUMBER;
	for(i=0;i<m;i++)
		{ts=gonl_arra_govj(temp,k,idx[i],&dummy);
		if(ts==0)
			{dis=petn_dist_vibr(C[s],C[idx[i]]);
			if(dis<sml)
				{sml=dis;
				q=i;
				}
			}
		}
	temp[k]=idx[q];
	}

for(i=0;i<m;i++)
	idx[i]=temp[i];
free(temp);
}



int rofq_extr_hopm(ns_curv *C,int n,
int *excl,int p,double eps,int *idx)
{int m,i,j,ts,tx,tr,dummy,m_new=0;

for(i=0;i<n;i++)
	{ts=gonl_arra_govj(excl,p,i,&dummy);
	if(ts==0)
		idx[0]=i;
	}
m=1;

for(j=0;j<n;j++)
	{for(i=0;i<n;i++)
		{tx=gonl_arra_govj(excl,p,i,&dummy);
		if(tx==0)
			{ts=wuhj_conn_bonk(C,idx,m,C[i],eps);
			tr=gonl_arra_govj(idx,m,i,&dummy);
			if((ts==1)&&(tr==0))
				{idx[m]=i;
				m++;
				}
			}
		}
	if(m_new==m)
		break;
	m_new=m;
	}
rijn_rear_woqp(C,n,idx,eps,m);
return m;
}


void fuch_find_tams(ns_curv *C,int n,double eps,
trm_ids *I)
{int *idx,*excl,p,m,i,j,k;
idx=(int *)malloc((n*n+20)*sizeof(int));
excl=(int *)malloc((n*n+20)*sizeof(int));
p=0;	k=0;
for(i=0;i<n;i++)
	{m=rofq_extr_hopm(C,n,excl,p,eps,idx);
	
	I->nb_curve_comp[k]=m;
	for(j=0;j<m;j++)
		{if(idx[j]<0)
			{fprintf(tmpout,"ng idx\n");
			exit(0);
			}
		I->list_curve_comp[k][j]=idx[j];
		}
	
	for(j=0;j<m;j++)
		excl[p+j]=idx[j];
	p=p+m;
	k++;
	if(p==n)
		break;
	}
free(idx);
free(excl);
I->nb_inter=k-1;
}



void wiqs_conn_vomw(ns_curv *C,trm_ids I,
c_curve *cc)
{int i,p,N,s;
for(p=0;p<I.nb_inter+1;p++)
	{N=I.nb_curve_comp[p];
	cc[p].N=N;
	cc[p].nle=0;
	cc[p].nca=0;
	cc[p].nnc=N;
	for(i=0;i<N;i++)
		{s=I.list_curve_comp[p][i];
		zobm_find_wumq_kihf(C[s],&cc[p].nc[i]);
		cc[p].type[i]=2;
		}
	quhw_flip_zoph(&cc[p]);
	}
}



int wong_comp_golp(trm_sph T,
c_arc3D *C_loc,int n,c_curve *cc)
{int m,i,n_loc;
double eps=0.01;
prop_n_curv pnc;
ns_curv *nc;
trm_ids I;

n_loc=n;
I.nb_curve_comp=(int *)malloc(n_loc*sizeof(int));
I.list_curve_comp=(int **)malloc(n_loc*sizeof(int*));
for(i=0;i<n_loc;i++)
	I.list_curve_comp[i]=(int *)malloc(n_loc*sizeof(int));

pnc.n=4;
pnc.k=3;
nc=(ns_curv *)malloc(n*sizeof(ns_curv));
for(i=0;i<n;i++)
	{foks_allo_vukp(pnc,&nc[i]);
	dolj_curv_kacq(T,C_loc[i],&nc[i]);
	}

fuch_find_tams(nc,n,eps,&I);
m=I.nb_inter;


wiqs_conn_vomw(nc,I,cc);
for(i=0;i<n;i++)
	newt_dest_lefq(pnc,&nc[i]);
free(nc);


for(i=0;i<n_loc;i++)
	free(I.list_curve_comp[i]);
free(I.list_curve_comp);
free(I.nb_curve_comp);
return m+1;
}

 
void jotc_find_jewl_zorv(int nb,prop_ccurve *pcc)
{int i;
pcc->nca=0;
pcc->nle=0;
pcc->nnc=nb;
pcc->N=nb;
for(i=0;i<nb;i++)
	{pcc->n[i]=4;
	pcc->k[i]=3;
	}
}


double medn_fart_rakq(parm *X,int N)
{int i;
double lrg,dis,x,y;
lrg=0.0;
for(i=0;i<N;i++)
	{x=X[i].u;
	y=X[i].v;
	dis=x*x+y*y;
	if(dis>lrg)
		lrg=dis;
	}
return lrg;
}


int tulj_curr_winq(polygon *P,int N,int *excl,int m)
{int q=-1,i,dummy,ts;
double lrg,sz;
lrg=0.0;
for(i=0;i<N;i++)
	{ts=gonl_arra_govj(excl,m,i,&dummy);
	if(ts==0)
		{sz=medn_fart_rakq(P[i].vertex,P[i].v_grs);
		if(sz>lrg)
			{lrg=sz;
			q=i;
			}
		}
	}
return q;
}



void pugj_crea_sotm(int r,int *list,int m,trm_sph T,
c_curve *cc,int n,trmsrf *S,int *orient)
{int ort,i,k;
S->nb_inner=n-1;
if(m>=MAX_INTERNAL_CURVES)
	{fprintf(tmpout,"MAX_INTERNAL_CURVES is reached\n");
	exit(0);
	}

if(cc[r].N>=MAXCOMP)
	{fprintf(tmpout,"2.  MAXCOMP=%d is reached\n",MAXCOMP);
	exit(0);
	}
ort=sufc_orie_cuvf(cc[r]);
if(ort==COUNTER_CLOCKWISE)
	{kotg_find_wuhk_kemt(cc[r],&S->cc);
	orient[r]=+1;
	}
if(ort==CLOCKWISE)
	{mevj_inve_nujf(cc[r],&S->cc);
	orient[r]=-1;
	}

S->nb_inner=m;
for(k=0;k<m;k++)
	{i=list[k];
	if(cc[i].N>=MAXCOMP)
		{fprintf(tmpout,"3.   MAXCOMP=%d is reached\n",MAXCOMP);
		exit(0);
		}
	ort=sufc_orie_cuvf(cc[i]);
	if(ort==CLOCKWISE)
		{kotg_find_wuhk_kemt(cc[i],&S->inner[k]);
		orient[i]=+1;
		}
	if(ort==COUNTER_CLOCKWISE)
		{mevj_inve_nujf(cc[i],&S->inner[k]);
		orient[i]=-1;
		}
	}

S->type=3;
zikt_find_jotz_jewb(T,&S->ts);
S->boundary=1;
}


int hect_list_mudc(sphere *S,int nb_sph,adj_hash H,
set_arcs *SA,trmsrf *surf,int *supp,int max_surf)
{int i,j,nb_cur,p,nb_comp,nb_arcs;
trm_sph *T;
prop_ccurve pcc;
c_curve *cc;
T=(trm_sph *)malloc(nb_sph*sizeof(trm_sph));
sonv_trim_gopl(S,nb_sph,T,H);
nb_cur=0;
for(i=0;i<nb_sph;i++)
if(SA[i].ar_grs!=1)
	{nb_arcs=SA[i].ar_grs;
	if(nb_arcs==1)
		{fprintf(tmpout,"Warning: one arc only on the atom %d/%d\n",i,nb_arcs-1);
		exit(0);
		}
	if(verbose_variable==VERBOSE)
		{if(((i%10)==0)||(i==nb_sph-1))
			fprintf(tmpout,"Trimming atom surf[%d/%d]   nb_cur=%d  nb_arcs=%d\n",i,nb_sph-1,nb_cur,nb_arcs);
		}
	cc=(c_curve *)malloc(nb_arcs*sizeof(c_curve));
	jotc_find_jewl_zorv(nb_arcs,&pcc);
	for(j=0;j<nb_arcs;j++)
		homd_allo_tevf(pcc,&cc[j]);
	nb_comp=wong_comp_golp(T[i],SA[i].C,nb_arcs,cc);
	p=qonk_form_sevt(i,T[i],nb_arcs,cc,nb_comp,
	surf,nb_cur,supp,max_surf);
	nb_cur=nb_cur+p;
	for(j=0;j<nb_arcs;j++)
		wosn_dest_jomw(pcc,&cc[j]);
	free(cc);
	}
free(T);
if(verbose_variable==VERBOSE)
	fprintf(tmpout,"trim search is complete  nb_cur=%d\n",nb_cur);
return nb_cur;
}



