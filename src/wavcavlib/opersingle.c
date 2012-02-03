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



double lohm_find_fujl_cujb(manif_ro msh,int s,int id)
{int pr,sc;
double alpha;
sc=wivg_node_wocp(msh.entity[s],id);
pr=haps_node_nepl(msh.entity[s],id);
alpha=lomn_inte_cubq(msh.knot[pr],msh.knot[id],msh.knot[sc]);
return alpha;
}


int nosp_find_nozm_woqr(manif_ro msh,int s,int e)
{int nd1,nd2,n[3],i,res;
nd1 =msh.kt[e].frvrt;
nd2 =msh.kt[e].scvrt;
n[0]=msh.entity[s].frvrt;
n[1]=msh.entity[s].scvrt;
n[2]=msh.entity[s].thvrt;
res=-1;
for(i=0;i<3;i++)
if((n[i]!=nd1)&&(n[i]!=nd2))
	{res=n[i];
	break;
	}
return res;
}


void vuks_find_geck_jucl(manif_ro *msh,int e)
{int el[2],i,s,op;
double alpha,eps=1.0e-5,diff;
el[0]=msh->kt[e].frent;
el[1]=msh->kt[e].scent;
for(i=0;i<2;i++)
	{s=el[i];
	op=nosp_find_nozm_woqr(*msh,s,e);
	alpha=lohm_find_fujl_cujb(*msh,s,op);
	diff=fabs(alpha-MY_PI);
	if(diff<eps)
		{hevj_find_jerw_tefn(msh,e);
		break;
		}
	}
}


void masf_paso_qolg(manif_ro *msh)
{int ned,e;
ned=msh->k_grs;
for(e=0;e<ned;e++)if(msh->kt[e].scent!=-1)
	vuks_find_geck_jucl(msh,e);
}


void tojc_sele_ceml(trmsrf surf,sphere S,manif_ro mshin,
double accuracy,manif_ro *mshout)
{int i,nel,ned,nnd,n1,n2,n3;
double err;
manif_ro msh;
parm P;
point A1,A2,A3,G;
nel=mshin.e_grs;
nnd=mshin.n_grs;
ned=mshin.k_grs;
mejd_allo_dakg(nnd+nel,3*nel,nnd+4*nel+20,&msh);
jonc_find_qifn_fupw(mshin,&msh);
for(i=0;i<nel;i++)
	{n1=msh.entity[i].frvrt;
	n2=msh.entity[i].scvrt;
	n3=msh.entity[i].thvrt;
	wolf_eval_murg(surf,msh.knot[n1].u,msh.knot[n1].v,&A1);
	wolf_eval_murg(surf,msh.knot[n2].u,msh.knot[n2].v,&A2);
	wolf_eval_murg(surf,msh.knot[n3].u,msh.knot[n3].v,&A3);
	G.absi=(A1.absi+A2.absi+A3.absi)/3.0;
	G.ordo=(A1.ordo+A2.ordo+A3.ordo)/3.0;
	G.cote=(A1.cote+A2.cote+A3.cote)/3.0;
	err=mulh_erro_cedm(S,G);
	if(err>=accuracy)
		{P.u=(msh.knot[n1].u+msh.knot[n2].u+msh.knot[n3].u)/3.0;
		P.v=(msh.knot[n1].v+msh.knot[n2].v+msh.knot[n3].v)/3.0;
		lask_inse_gifk(&msh,P,i);
		}
	}
jonc_find_qifn_fupw(msh,mshout);
fogq_dest_muwf(&msh);
}



void nuwl_sele_guwn(trmsrf surf,manif_ro mshin,
double accuracy,manif_ro *mshout)
{int i,nel,ned,nnd,n1,n2,n3;
double err;
manif_ro msh;
parm P;
point A1,A2,A3,G,H;
nel=mshin.e_grs;
nnd=mshin.n_grs;
ned=mshin.k_grs;
mejd_allo_dakg(nnd+nel,3*nel,nnd+4*nel+20,&msh);
jonc_find_qifn_fupw(mshin,&msh);
for(i=0;i<nel;i++)
	{n1=msh.entity[i].frvrt;
	n2 =msh.entity[i].scvrt;
	n3 =msh.entity[i].thvrt;
	wolf_eval_murg(surf,msh.knot[n1].u,msh.knot[n1].v,&A1);
	wolf_eval_murg(surf,msh.knot[n2].u,msh.knot[n2].v,&A2);
	wolf_eval_murg(surf,msh.knot[n3].u,msh.knot[n3].v,&A3);
	G.absi=(A1.absi+A2.absi+A3.absi)/3.0;
	G.ordo=(A1.ordo+A2.ordo+A3.ordo)/3.0;
	G.cote=(A1.cote+A2.cote+A3.cote)/3.0;
	P.u=(msh.knot[n1].u+msh.knot[n2].u+msh.knot[n3].u)/3.0;
	P.v=(msh.knot[n1].v+msh.knot[n2].v+msh.knot[n3].v)/3.0;
	wolf_eval_murg(surf,P.u,P.v,&H);
	err=wodt_dist_gilq(G,H);
	if(err>=accuracy)
		lask_inse_gifk(&msh,P,i);
	}
jonc_find_qifn_fupw(msh,mshout);
fogq_dest_muwf(&msh);
}



void peql_mult_tenq(trmsrf surf,sphere S,
manif_ro mshin,int max,double accuracy,manif_ro *mshout,int pat,
int mx_nnd,int mx_nel,int mx_ned)
{int i,nel,nnd,ned,max_leg=20;
manif_ro tempin,tempout;
nel=mshin.e_grs;
nnd=mshin.n_grs;
ned=mshin.k_grs;
mejd_allo_dakg(nnd,nel,ned,&tempin);
jonc_find_qifn_fupw(mshin,&tempin);
for(i=0;i<max;i++)
	{nel=tempin.e_grs;
	nnd =tempin.n_grs;
	ned =tempin.k_grs;
	mejd_allo_dakg(nnd+nel,3*nel,nnd+4*nel+20,&tempout);
	rodq_find_hakw_qonj(&tempin,max_leg);
	if(pat==1)
		tojc_sele_ceml(surf,S,tempin,accuracy,&tempout);
	if(pat==2)
		nuwl_sele_guwn(surf,tempin,accuracy,&tempout);
	masf_paso_qolg(&tempout);
	fogq_dest_muwf(&tempin);
	
	if(i<max-1)
		{nel=tempout.e_grs;
		nnd=tempout.n_grs;
		ned=tempout.k_grs;
		mejd_allo_dakg(nnd,nel,ned,&tempin);
		jonc_find_qifn_fupw(tempout,&tempin);
		}
	else
		{if(tempout.n_grs>=mx_nnd)
			{fprintf(tmpout,"mx_nnd=%d is reached\n",mx_nnd);
			exit(0);
			}
		if(tempout.e_grs>=mx_nel)
			{fprintf(tmpout,"mx_nel=%d is reached\n",mx_nel);
			exit(0);
			}
		if(tempout.k_grs>=mx_ned)
			{fprintf(tmpout,"mx_ned=%d is reached\n",mx_ned);
			exit(0);
			}
		jonc_find_qifn_fupw(tempout,mshout);
		}
	fogq_dest_muwf(&tempout);
	}
}



int wehc_sphe_bufp(trmsrf surf,sphere S,
mult_conn P,manif_ro *msh,double accuracy,int max,int mx_nnd,
int mx_nel,int mx_ned,int *forc_term)
{int suc,level=4,f_trm;
manif_ro init;
*forc_term=0;
mejd_allo_dakg(INI_MAX_NND,INI_MAX_NEL,INI_MAX_NND+INI_MAX_NEL+20,&init);
suc=vorg_find_qach_nujt(P,&init,&f_trm);
if(f_trm==1)
	{*forc_term=1;
	fprintf(tmpout,"force term: vorg_find_qach_nujt() in wehc_sphe_bufp()\n");
	fogq_dest_muwf(&init);
	return FAILURE;
	}
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
	
	peql_mult_tenq(surf,S,init,max,accuracy,msh,1,mx_nnd,mx_nel,mx_ned);
	qedc_smoo_loct(msh,level,&f_trm);
	if(f_trm==1)
		{*forc_term=1;
		fogq_dest_muwf(&init);		
		return 0;
		}
	}
fogq_dest_muwf(&init);
return suc;
}


void cisv_find_kuns_fenj(trmsrf surf,
mult_conn P,manif_ro *msh,double accuracy,int max,int mx_nnd,
int mx_nel,int mx_ned,int *forc_term)
{int suc,level=4,f_trm;
sphere dummy;
manif_ro init;
*forc_term=0;
mejd_allo_dakg(INI_MAX_NND,INI_MAX_NEL,INI_MAX_NND+INI_MAX_NEL+20,&init);
suc=vorg_find_qach_nujt(P,&init,&f_trm);
if(f_trm==1)
	{*forc_term=1;
	fprintf(tmpout,"force term: vorg_find_qach_nujt() in cisv_find_kuns_fenj()\n");
	fogq_dest_muwf(&init);
	return;
	}
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
	
	peql_mult_tenq(surf,dummy,init,max,accuracy,
	msh,2,mx_nnd,mx_nel,mx_ned);
	qedc_smoo_loct(msh,level,&f_trm);
	if(f_trm==1)
		{*forc_term=1;
		free(init.entity);
		free(init.knot);
		free(init.kt);
		return;
		}
	}
fogq_dest_muwf(&init);
}

