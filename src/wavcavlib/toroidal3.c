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
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"
#include "sas.h"


int fiws_inte_puds(parent_circle *PC,adj_hash H,sphere *S,
int nb_sph,double probe,int w1,int w2,set_arcs *SA,pt_tor *PT,
int *arc_id1,int *arc_id2,int *map1,int *map2,int max_pos)
{int nb1,nb2,*ar1,*ar2,m1,m2;
int n_loc,ts;

m1=SA[w1].ar_grs;
m2=SA[w2].ar_grs;
ar1=(int *)malloc(m1*sizeof(int));
ar2=(int *)malloc(m2*sizeof(int));
ts=mefj_list_hepq(PC,H,w1,w2,S,nb_sph,probe,SA,ar1,ar2,&nb1,&nb2);
if(ts==0)
	{free(ar1);
	free(ar2);
	return 0;
	}


if(nb1!=nb2)
	{fprintf(tmpout,"expected?  [%d,%d]\n",nb1,nb2);
	}
n_loc=gedh_poss_jucp(S,probe,w1,w2,SA,ar1,ar2,
nb1,nb2,PT,arc_id1,arc_id2,max_pos);

free(ar1);		free(ar2);
return n_loc;
}


void rakg_veri_cizm(sphere *S,pt_tor *PT,int nb,
set_arcs *SA,int w1,int w2,double eps)
{int j;
double d1,d2;
for(j=0;j<nb;j++)
	{d1=vuqg_dist_faql(PT[j],S[w1]);
	d2=vuqg_dist_faql(PT[j],S[w2]);
	if((d1>eps)||(d2>eps))
		{fprintf(tmpout,"Warning: mismatch\n");
		exit(0);
		}
	}
fprintf(tmpout,"Good local match\n");
}



int luhb_inco_zokt(parent_circle *PC,adj_hash H,
double probe,sphere *S,int nb_sph,set_arcs *SA,
int nb_cur,trmsrf *surf,blend_cpx *BC,int max_surf)
{int w1,w2,nb,nb_new,j,k,m,dummy_i,dummy_j,ts,q;
int newpt,*arc_id1,*arc_id2,*map1,*map2;
double eps_paral=1.0e-4;
pt_tor *PT;
k=nb_cur;	newpt=0;
for(w1=0;w1<nb_sph;w1++)
	{if(verbose_variable==VERBOSE)
		fprintf(tmpout,"Smoothing with atom: %d / %d \n",w1,nb_sph-1);
	for(w2=0;w2<w1;w2++)
		{ts=qoml_find_venq_fugh(probe,S,SA,w1,w2,&dummy_i,&dummy_j,eps_paral);
		if(ts==0)
			{m=SA[w1].ar_grs+SA[w2].ar_grs+10;
			PT=(pt_tor *)malloc(m*sizeof(pt_tor));
			arc_id1=(int *)malloc(m*sizeof(int));
			arc_id2=(int *)malloc(m*sizeof(int));
			map1=(int *)malloc(m*sizeof(int));
			map2=(int *)malloc(m*sizeof(int));
			
			nb=fiws_inte_puds(PC,H,S,nb_sph,probe,w1,w2,SA,
			PT,arc_id1,arc_id2,map1,map2,m);
			if(nb>=m) 
				{fprintf(tmpout,"maximal allocation reached\n");
				exit(0);
				}
			
			
			
			for(j=0;j<nb;j++)
				{surf[k].type=4;
				surf[k].boundary=1;
				romh_find_cont_qucr(PT[j],&surf[k].pt);
				sart_rect_jamc(0.0,1.0,0.0,1.0,&surf[k].cc);
				surf[k].nb_inner=0;
				
				newpt++;
				
				q=BC->bt_grs;
				if(q>=MAX_BLEND_TOR)
					{fprintf(tmpout,"MAX_BLEND_TOR is reached\n");
					exit(0);
					}
				BC->BT[q].sph_idx1=w1;
				BC->BT[q].sph_idx2=w2;
				BC->BT[q].trim_idx=k;
				BC->BT[q].a_id1=arc_id1[j];
				BC->BT[q].a_id2=arc_id2[j];
				BC->bt_grs=q+1;
				
				k++;
				
				if(k>=max_surf)
					{fprintf(tmpout,"Max nb trimmed surfaces is reached\n");
					exit(0);
					}
				}
			free(PT);
			free(map1);
			free(map2);
			free(arc_id1);
			free(arc_id2);
			}
		}
	}
nb_new=k;
fprintf(tmpout,"Current nb trimmed surf=%d\n",nb_new);
fprintf(tmpout,"New patches=%d\n",newpt);

return nb_new;
}


 
