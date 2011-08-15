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
#include <malloc.h>
#include <stdlib.h>
#include <math.h> 
#include "cavity.h"  
#include "pln_sph.h"
#include "sas.h"
#include "splinemol.h"
#include "meshsas.h"
#include "partcas.h"
#include "geodesic.h"
#include "smooth.h"


double	prec_probe;
int		prec_ex;



void bugt_clea_pujv(supp_sph_tri *treb,int nb_surf3)
{int i;
supp_sph_tri *treb_temp;
treb_temp=(supp_sph_tri *)malloc(nb_surf3*sizeof(supp_sph_tri));
for(i=0;i<nb_surf3;i++)
	fewd_find_lepz_zaqn(treb[i],&treb_temp[i]);
free(treb);
treb=(supp_sph_tri *)malloc(nb_surf3*sizeof(supp_sph_tri));
for(i=0;i<nb_surf3;i++)
	fewd_find_lepz_zaqn(treb_temp[i],&treb[i]);
free(treb_temp);
}


void tofr_sort_mebd(atom *A,adj_hash H,int p)
{int N1,N2,i,w;
double x,y,z,r;
FILE *fp;
fprintf(tmpout,"Atom[%d]\n",p);
N1=H.entry[p].nb_neighbors;
N2=H.inter[p].nb_neighbors;
fprintf(tmpout,"size:  adjacency=%d   interaction=%d\n",N1,N2);
fp=fopen("sort.dat","w");
x=A[p].zent.absi;
y=A[p].zent.ordo;
z=A[p].zent.cote;
r=A[p].rad;
fprintf(fp,"%.20f   %.20f   %.20f   %.20f\n",x,y,z,r);
for(i=0;i<N1;i++)
	{w=H.entry[p].neighbor[i];
	x=A[w].zent.absi;
	y=A[w].zent.ordo;
	z=A[w].zent.cote;
	r=A[w].rad;
	fprintf(fp,"%.20f   %.20f   %.20f   %.20f\n",x,y,z,r);
	}
for(i=0;i<N2;i++)
	{w=H.inter[p].neighbor[i];
	x=A[w].zent.absi;
	y=A[w].zent.ordo;
	z=A[w].zent.cote;
	r=A[w].rad;
	fprintf(fp,"%.20f   %.20f   %.20f   %.20f\n",x,y,z,r);
	}
fclose(fp);
}


void hetc_glob_kaqc(double probe,double c_param,sphere *S,
int nb_sph,adj_hash H,int *forc_term)
{int i,ans,nb_pc,new_pc,max_lev_mis=1,*V1,*V2,nb_V;
int max_tor=MAX_TOR_PATCH,max_self=MAX_SELFINTERS;
int nb_surf1=0,nb_surf2=0,nb_surf3=0,new_surf,*supp,meth=3;
int n_gap,max_gap=100,max_gap_val=20,f_trm,f_trm2,f_trm3;
double thres_stretch=0.1;
trmsrf *surf1,*surf2,*surf3;
supp_sph_tri *treb;
parent_circle *PC;
blend_cpx BC;
hash_entry *gap;
set_arcs *SA;
*forc_term=0;
fprintf(tmpout,"Organizing B-Rep\n");
PC=(parent_circle *)malloc(nb_sph*nb_sph*sizeof(parent_circle));
SA=(set_arcs *)malloc(nb_sph*sizeof(set_arcs));
nb_pc=0;
for(i=0;i<nb_sph;i++)
	{if(verbose_variable==VERBOSE)
		fprintf(tmpout,"Processing atom [%d/%d]\n",i,nb_sph-1);
	SA[i].C=(c_arc3D *)malloc(MAX_ARCS*sizeof(c_arc3D));
	SA[i].par_idx=(int *)malloc(MAX_ARCS*sizeof(int));
	new_pc=hagf_arcs_belz(H,probe,S,i,&SA[i],PC,nb_pc);
	nb_pc=new_pc;
	}
gojw_pair_nokv(nb_sph,SA,PC,nb_pc);
raqn_inve_kotj(nb_sph,SA);
surf1=(trmsrf *)malloc(MAX_TRIM_SURF1*sizeof(trmsrf));
hegp_allo_qogc(surf1,MAX_TRIM_SURF1,MAX_INTERNAL_CURVES);
supp=(int *)malloc(MAX_TRIM_SURF1*sizeof(int));
fprintf(tmpout,"Forming trimmed surfaces\n");
nb_surf1=hect_list_mudc(S,nb_sph,H,SA,surf1,supp,MAX_TRIM_SURF1);
raqn_inve_kotj(nb_sph,SA);
zaqw_dest_jelq(surf1,nb_surf1,MAX_TRIM_SURF1,MAX_INTERNAL_CURVES);
surf2=(trmsrf *)malloc(MAX_TRIM_SURF2*sizeof(trmsrf));
culj_allo_pudn(surf2,MAX_TRIM_SURF2);
V1=(int *)malloc(max_tor*sizeof(int));
V2=(int *)malloc(max_tor*sizeof(int));
fprintf(tmpout,"Find some toroidal\n");
nb_surf2=sapn_some_judq(meth,surf1,nb_surf1,supp,PC,&nb_pc,H,
S,nb_sph,probe,SA,surf2,0,max_tor,max_self,MAX_TRIM_SURF2,
V1,V2,thres_stretch,&f_trm);
if(f_trm==1)
	{fprintf(tmpout,"force term: sapn_some_judq() in hetc_glob_kaqc()\n");
	*forc_term=1;
	}
if(f_trm==0)
	{mehz_chec_dijp(nb_sph,SA);
	nb_V=nb_surf2;
	if(verbose_variable==VERBOSE)
		fprintf(tmpout,"Current number of trimmed surfs=%d\n",nb_surf1+nb_surf2);
	raqn_inve_kotj(nb_sph,SA);
	qatn_allo_mobc(nb_sph,&BC);
	fprintf(tmpout,"Find incomplete toroidal\n");
	new_surf=luhb_inco_zokt(PC,H,probe,S,nb_sph,
	SA,nb_surf2,surf2,&BC,MAX_TRIM_SURF2);
	raqn_inve_kotj(nb_sph,SA);
	}
free(PC);
if(f_trm==0)
	{nb_surf2=new_surf;
	if(verbose_variable==VERBOSE)
		{fprintf(tmpout,"nb_surf2=%d\n",nb_surf2);
		fprintf(tmpout,"Current number of trimmed surfs=%d\n",nb_surf1+nb_surf2);
		}
	kejh_fill_zogq(S,surf2,nb_sph,&BC);
	fprintf(tmpout,"Update merge:\n");
	if(verbose_variable==VERBOSE)
		fprintf(tmpout,"Fill patch concave\n");
	surf3=(trmsrf *)malloc(MAX_TRIM_SURF3*sizeof(trmsrf));
	gilp_allo_temc(surf3,MAX_TRIM_SURF3);
	treb=(supp_sph_tri *)malloc(MAX_TRIM_SURF3*sizeof(trmsrf));
	nb_surf3=hicj_fill_zalj(probe,S,nb_sph,surf2,0,BC,surf3,treb,MAX_TRIM_SURF3);
	legw_upda_tevg(nb_sph,&BC,V1,V2,nb_V);
	for(i=0;i<BC.bt_grs;i++)
		vupn_chec_wudq(surf2,S,BC.BT[i]);
	}
free(V1);	free(V2);
if(f_trm==0)
	{fprintf(tmpout,"nb_blend_tor=%d\n",BC.bt_grs);
	new_surf=hefk_apen_sikf(probe,6,S,nb_sph,surf1,nb_surf1,supp,
	surf2,nb_surf2,&BC,SA,MAX_TRIM_SURF2);
	nb_surf2=new_surf;
	novw_dest_wojm(surf2,nb_surf2,MAX_TRIM_SURF2);
	fprintf(tmpout,"Spherical representation\n");
	homg_sphe_vemp(probe,surf3,nb_surf3);
	fprintf(tmpout,"FIX SELF-INTERSECTION\n");
	fprintf(tmpout,"Current number of trimmed surfs=%d\n",nb_surf1+nb_surf2+nb_surf3);
	fprintf(tmpout,"BOUNDS: [%d,%d,%d]\n",nb_surf1,nb_surf2,nb_surf3);
	kejh_fill_zogq(S,surf2,nb_sph,&BC);
	zevt_stru_copb(H,S,nb_sph,SA,supp,surf1,nb_surf1,surf2,BC,&f_trm3);
	if(f_trm3==1)
		{*forc_term=1;
		fprintf(tmpout,"force term: zevt_stru_copb() in hetc_glob_kaqc()\n");
		}
	if(f_trm3==0)
		{gap=(hash_entry *)malloc(max_gap*sizeof(hash_entry));
		for(i=0;i<max_gap;i++)
			gap[i].list=(int *)malloc(max_gap_val*sizeof(int));
		new_surf=qumw_atom_picm(S,probe,surf1,nb_surf1,surf2,nb_surf2,
		surf3,nb_surf3,SA,BC,1.0e-4,MAX_TRIM_SURF3,gap,&n_gap,max_gap,max_gap_val);
		nb_surf3=new_surf;
		}
	luqw_dest_horb(surf3,nb_surf3,MAX_TRIM_SURF3);
	if(f_trm3==0)
		{fprintf(tmpout,"n_gap=%d\n",n_gap);
		if(verbose_variable==VERBOSE)
			{fprintf(tmpout,"Current number of trimmed surfs=%d\n",nb_surf1+nb_surf2+nb_surf3);
			  /*
			    fprintf(tmpout,"Dump molecular surface? ans=");
			    scanf("%d",&ans);
			    if(ans==1)
			  */
			  dump_molc_surf(surf1,nb_surf1,surf2,nb_surf2,surf3,nb_surf3);
			}
		luvm_blen_foms(S,nb_sph,supp,surf1,nb_surf1,
		surf2,nb_surf2,surf3,nb_surf3,SA,BC,gap,n_gap,c_param,&f_trm2);
		if(f_trm2==1)
			{fprintf(tmpout,"force term: luvm_blen_foms() in hetc_glob_kaqc()\n");
			*forc_term=1;
			}
		for(i=0;i<max_gap;i++)
			free(gap[i].list);
		free(gap);
		luqw_dest_horb(surf3,0,nb_surf3);
		}
	free(surf3);
	free(treb);
	bulf_dest_tacs(nb_sph,&BC);
	}
for(i=0;i<nb_sph;i++)
	{free(SA[i].C);
	free(SA[i].par_idx);
	}
free(SA);
zaqw_dest_jelq(surf1,0,nb_surf1,MAX_INTERNAL_CURVES);
free(surf1);
novw_dest_wojm(surf2,0,nb_surf2);
free(surf2);
free(supp);
}


int bihp_test_romg(char *filename, double probe,
double c_param,int *n_size)
{int nb_sph,new_nb,*list,res;
int i,z,ts,dummy,nb_I,k,f_trm;
int maxval=100,max_inter=400;
atom *S,*I,*R;
adj_hash H;
if(fabs(probe)<1.0e-7)
	{fprintf(tmpout,"Zero probe radius\n");
	exit(0);
	}
nb_sph=vofr_numb_qimf(filename,&S); // Number of spheres (from "ex")
if(nb_sph < 0) return -1;
*n_size=nb_sph;
fprintf(tmpout,"SIZE=%d atoms\n",nb_sph);
fprintf(tmpout,"Number of spheres=%d\n",nb_sph);
list=(int *)malloc(nb_sph*sizeof(int));
fprintf(tmpout,"Isolated\n");
nb_I=varm_dete_qumz(S,nb_sph,probe,list);
if(nb_I>=1)
	{I=(atom *)malloc(nb_I*sizeof(atom));
	for(i=0;i<nb_I;i++)
		{z=list[i];
		neqg_find_lodr_bogm(S[z],&I[i]);
		}
	sotp_atom_kegt(I,nb_I);
	free(I);
	}
R=(atom *)malloc((nb_sph-nb_I)*sizeof(atom));
k=0;
for(i=0;i<nb_sph;i++)
	{ts=gonl_arra_govj(list,nb_I,i,&dummy);
	if(ts==0)
		{neqg_find_lodr_bogm(S[i],&R[k]);
		k++;
		}
	}
free(list);
nb_sph=k;
free(S);
fprintf(tmpout,"Removal of nestedness\n");
new_nb=zosd_remo_qohj(R,nb_sph);
nb_sph=new_nb;
new_nb=qefl_remo_gopd(R,nb_sph);
nb_sph=new_nb;
if(nb_sph>MAX_ATOMS)
	{fprintf(tmpout,"Number of atoms=%d\n",nb_sph);
	fprintf(tmpout,"Maximum size of molecule=%d\n",MAX_ATOMS);
	exit(0);
	}
gikf_allo_nibl(nb_sph,maxval,max_inter,&H);
fprintf(tmpout,"FIND ADJACENCY:\n");
tojw_hash_pojh(probe,R,nb_sph,&H,maxval,max_inter);
hetc_glob_kaqc(probe,c_param,R,nb_sph,H,&f_trm); 
mewh_dest_selk(nb_sph,&H);
free(R); 
if(f_trm==0)
	{fprintf(tmpout,"GOOD TERMINATION\n");
	res=1;
	}
if(f_trm==1)
	{fprintf(tmpout,"PROBLEMS WERE DETECTED (FORCED QUIT)\n");
	res=0;
	}
return res;
}
 
extern int PATCH_LEVEL;

/*
int main(int argc,char **argv)
{
  int max_mod=33, dummy, test;
  double probe, c_param;
  verbose_variable=VERBOSE;
  FILE *tmpout = stdout;

  if(argc != 5){
    fprintf(stderr,"Usage: test.x INPUT_FILE.XYZR PROBE_RADIUS COARSITY_PARAMETER PATCH_LEVEL\n\n");
    fprintf(stderr,"Notes: Radii are not multiplied by 1.2.\n");
    fprintf(stderr,"       Coarsity is from 0.0 (coarse) to 1.0 (fine).\n\n");
    return -1;
  }
  
  char *endptr;
  probe = strtod(argv[2],&endptr);
  if(endptr ==  argv[2]){
    fprintf(stderr,"PROBE_RADIUS not good!\n");
    return -2;
  }

  c_param = strtod(argv[3],&endptr);
  if(endptr ==  argv[3]){
    fprintf(stderr,"COARSITY_PARAMETER not good!\n");
    return -3;
  }

  PATCH_LEVEL = strtol(argv[4],&endptr,10);
  if(endptr ==  argv[4]){
    fprintf(stderr,"PATCH_LEVEL not good!\n");
    return -3;
  }

  freopen("test.out","w",stdout);
  test = bihp_test_romg(argv[1],probe,c_param,&dummy);
  fclose(stdout);
  stdout = tmpout;

  if(test == 1){
    return 0; // GOOD
  }else{
    return -34; // BAD
  }
}
*/

FILE *tmpout;

 /*
  * @brief Interface routine to the cavity generator
  *
  * NOTE: Radii in "cavity.inp" are not multiplied by 1.2
  *
  * @param probe Probe radius
  * @param coarsity Coarsity parameter ]0.0, 1.0[
  * @param pl Patch level
  * @param info Changed to 0 if all went fine.
  *
  */
void cavity_create_(double *probe, double *coarsity, int *pl, int *info){
  char *infile = "cavity.inp";
  int dummy, test;

  if(*probe < 0.0){
    *info = -1;
    return;
  }
  if(*coarsity <= 0.0 || *coarsity >= 1.0){
    *info = -2;
    return;
  }
  if(*pl < 0){
    *info = -3;
    return;
  }

  PATCH_LEVEL = *pl;
  verbose_variable=VERBOSE;


  tmpout = fopen("create_cavity.out","w");
  test = bihp_test_romg(infile,*probe,*coarsity,&dummy);
  fclose(tmpout);

  if(test == 1){
    *info = 0;
  } else {
    *info = -4;
  }

  return;
};
  

