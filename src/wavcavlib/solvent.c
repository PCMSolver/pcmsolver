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
#include <malloc.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"
#include "sas.h"
#include "splinemol.h"
#include "meshsas.h"



void solr_find_vutk_dipl(manif_ro msh,trmsrf surf,manif_tl *MSH)
{int i,nnd,nel;
double x,y;
nnd=msh.n_grs;
nel=msh.e_grs;

for(i=0;i<nnd;i++)
	{x=msh.knot[i].u;
	y=msh.knot[i].v;
	wolf_eval_murg(surf,x,y,&MSH->knot[i]);
	}
for(i=0;i<nel;i++)
	{MSH->entity[i].frvrt=msh.entity[i].frvrt;
	MSH->entity[i].scvrt=msh.entity[i].scvrt;
	MSH->entity[i].thvrt=msh.entity[i].thvrt;
	}
MSH->n_grs=nnd;
MSH->e_grs=nel;
}


void mehz_chec_dijp(int nb_sph,set_arcs *SA)
{int i,j;
for(i=0;i<nb_sph;i++)
for(j=0;j<SA[i].ar_grs;j++)
if(SA[i].C[j].c_cir==1)
	{fprintf(tmpout,"Exist complet circle\n");
	exit(0);
	}
fprintf(tmpout,"No periodic CA\n");
}
 

void legw_upda_tevg(int nb_sph,blend_cpx *BC,
int *V1,int *V2,int nb_V)
{int N,i;
N=BC->bt_grs;
fprintf(tmpout,"Updating BC for sapn_some_judq  N=%d\n",N);
for(i=0;i<nb_V;i++)
	{if(N>=MAX_BLEND_TOR)
		{fprintf(tmpout,"MAX_BLEND_TOR is exceeded\n");
		exit(0);
		}
	BC->BT[N].a_id1=-1;
	BC->BT[N].a_id2=-1;
	BC->BT[N].sph_idx1=V1[i];
	BC->BT[N].sph_idx2=V2[i];
	BC->BT[N].trim_idx=i;
	N++;
	}
BC->bt_grs=N;
fprintf(tmpout,"Fill blend complex\n");
sedr_fill_luch(nb_sph,BC);
fprintf(tmpout,"Good blend complex\n");
}


void tebw_find_pokt_cifn(manif_ro msh,manif_ro *out)
{int nnd,nel,i;
nnd=msh.n_grs;
nel=msh.e_grs;
out->n_grs=nnd;
out->e_grs=nel;
for(i=0;i<nnd;i++)
	{out->knot[i].u=msh.knot[i].u;
	out->knot[i].v=msh.knot[i].v;
	}
for(i=0;i<nel;i++)
	{out->entity[i].frvrt=msh.entity[i].frvrt;
	out->entity[i].scvrt=msh.entity[i].scvrt;
	out->entity[i].thvrt=msh.entity[i].thvrt;
	}
}



void dahf_mega_zubj(manif_ro unit_quad,manif_ro unit_tri,
trmsrf *ts,int nb,int nb_all,megamanif *MG)
{int i,j,k,s,nin,depth=6,nb_awry,m[3],ind;
int tw,nnd=0,nel=0,N_quad=6,M_quad=6,f_trm;
double accuracy=0.003,offset_val=1.0e-4;
double *tm,mu_offset=0.5,eps_offset=0.5;
polygon P;
manif_ro msh;
parm omega;
sphere S;
mult_conn mc;

mejd_allo_dakg(IGS_MAX_NND,IGS_MAX_NEL,IGS_MAX_NED,&msh);
nb_awry=0;
k=MG->mw_grs;
for(j=0;j<nb;j++)
	{fprintf(tmpout,"1-TRMSURF=%d / %d    TP=%d\n",k,nb_all-1,ts[j].type);
	s=vubn_wors_qutk(ts[j].cc.N,depth);
	if((ts[j].type==4)||(ts[j].type==6))
		tebw_find_pokt_cifn(unit_quad,&msh);
	else if((ts[j].type==5)&&(ts[j].pc.sitn!=2))
		tebw_find_pokt_cifn(unit_tri,&msh);
	else if((ts[j].type==5)&&(ts[j].pc.sitn==2))
		{m[0]=10;	m[1]=10;	m[2]=10;
		welc_allo_dubg(0,m[0]+m[1]+m[2],&mc);
		gick_mult_jufs(m,ts[j].pc,&mc);
		getf_find_rogc_todj(ts[j].pc.T.zent,&S.zent);
		S.rad=ts[j].pc.T.rad;
		wehc_sphe_bufp(ts[j],S,mc,&msh,0.01,10,
		IGS_MAX_NND,IGS_MAX_NEL,IGS_MAX_NED,&f_trm);
		saqw_dest_kiqf(&mc);
		}
	else
		{nin=ts[j].nb_inner;
		ind=1;
		if(nin==0)
			{P.nb_local_vertices=(int *)malloc((nin+1)*sizeof(int));
			for(i=0;i<nin;i++) 
				s=s+vubn_wors_qutk(ts[j].inner[i].N,depth);
			P.vertex=(parm *)malloc(2*s*sizeof(parm));
			tm=(double *)malloc(2*s*sizeof(double));
			create_polygon(ts[j],&P,depth,accuracy,10,tm);
			tw=riwf_chec_zibc(P,1.0e-3,&omega);
			if(tw==1)
				{cojk_grad_negl(P.vertex,P.v_grs,omega,20,&msh);
				ind=2;
				}
			free(tm);
			free(P.vertex);
			free(P.nb_local_vertices);
			}
		if(ind==1)
			{welc_allo_dubg(nin,2*s,&mc);
			kanq_crea_tuqn(ts[j],&mc,depth,accuracy,10);
			getf_find_rogc_todj(ts[j].ts.zent,&S.zent);
			S.rad=ts[j].ts.rad;
			wehc_sphe_bufp(ts[j],S,mc,&msh,0.01,3,
			IGS_MAX_NND,IGS_MAX_NEL,IGS_MAX_NED,&f_trm);
			saqw_dest_kiqf(&mc);
			}
		if(msh.n_grs==0)
			{msh.e_grs=0;
			nb_awry++;
			}
		}
	MG->msh[k].knot=(point *)malloc(msh.n_grs*sizeof(point));
	MG->msh[k].entity=(telolf *)malloc(msh.e_grs*sizeof(telolf));
	solr_find_vutk_dipl(msh,ts[j],&MG->msh[k]);
	nnd=nnd+MG->msh[k].n_grs;
	nel=nel+MG->msh[k].e_grs;
	k++;
	}
fogq_dest_muwf(&msh);
MG->mw_grs=k;
fprintf(tmpout,"Number of awry cases=%d\n",nb_awry);
}


void dump_molc_surf(trmsrf *surf1,int nb1,
trmsrf *surf2,int nb2,trmsrf *surf3,int nb3)
{int j,nnd_tri,nel_tri;
int n_pat,p,N,M,i;
int nnd=0,nel=0,ans,N_quad=6,M_quad=6,N_tri=5;
manif_ro unit_quad,unit_tri;
megamanif MG;
FILE *fp;
unit_quad.knot=(parm *)malloc(N_quad*M_quad*sizeof(parm));
unit_quad.entity=(telolf *)malloc(2*(N_quad-1)*(M_quad-1)*sizeof(telolf));
lenq_simp_socg(N_quad,M_quad,&unit_quad);
nnd_tri=3*N_tri*N_tri;
nel_tri=6*(N_tri-1)*(N_tri-1);
unit_tri.knot=(parm *)malloc(nnd_tri*sizeof(parm));
unit_tri.entity=(telolf *)malloc(nel_tri*sizeof(telolf));
jatw_unit_hukl(N_tri,&unit_tri);
MG.msh=(manif_tl *)malloc((nb1+nb2+nb3)*sizeof(manif_tl));
MG.col=(rgb_lk *)malloc((nb1+nb2+nb3)*sizeof(rgb_lk));
MG.mw_grs=0;
dahf_mega_zubj(unit_quad,unit_tri,surf1,nb1,nb1+nb2+nb3,&MG);
dahf_mega_zubj(unit_quad,unit_tri,surf2,nb2,nb1+nb2+nb3,&MG);
dahf_mega_zubj(unit_quad,unit_tri,surf3,nb3,nb1+nb2+nb3,&MG);
free(unit_quad.knot);
free(unit_quad.entity);
free(unit_tri.knot);
free(unit_tri.entity);
n_pat=MG.mw_grs;
fp=fopen("dumped1.dat","w");
fprintf(fp,"%d\n",n_pat);
for(p=0;p<n_pat;p++)
  {N=MG.msh[p].n_grs;
  M=MG.msh[p].e_grs;
  fprintf(fp,"%d\n",N);
  fprintf(fp,"%d\n",M);
  for(i=0;i<N;i++)
    fprintf(fp,"%f  %f  %f\n",MG.msh[p].knot[i].absi,
    MG.msh[p].knot[i].ordo,MG.msh[p].knot[i].cote);
  for(i=0;i<M;i++)
    {fprintf(fp,"%d  ",MG.msh[p].entity[i].frvrt);
    fprintf(fp,"%d  ",MG.msh[p].entity[i].scvrt);
    fprintf(fp,"%d\n",MG.msh[p].entity[i].thvrt);
    }
  }
fclose(fp);
for(j=0;j<nb1+nb2+nb3;j++)
	{free(MG.msh[j].knot);
	free(MG.msh[j].entity);
	}
free(MG.msh);
free(MG.col);
fprintf(tmpout,"EXPORTED IN dumped1.dat\n");
}
