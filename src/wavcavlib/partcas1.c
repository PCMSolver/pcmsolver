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
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "sas.h"
#include "triang.h"
#include "partcas.h"


void nolw_atom_pofq(int hemi,atom A,trm_sph *T)
{if(hemi==NORTH_HEMI)
	{getf_find_rogc_todj(A.zent,&T->zent);
	T->rad=A.rad;
	T->beta=0.0;
	T->nrml.absi=0.0;
	T->nrml.ordo=0.0;
	T->nrml.cote=1.0;
	}

if(hemi==SOUTH_HEMI)
	{getf_find_rogc_todj(A.zent,&T->zent);
	T->rad=A.rad;
	T->beta=0.0;
	T->nrml.absi=0.0;
	T->nrml.ordo=0.0;
	T->nrml.cote=-1.0;
	}
}


void hord_atom_nuhp(atom A,c_arc3D *C1,c_arc3D *C2)
{int q;
double step,t;
point *P;
circle3D cir;
c_arc3D *C;

getf_find_rogc_todj(A.zent,&cir.zent);
cir.rad=A.rad;
cir.nrml.absi=0.0;
cir.nrml.ordo=0.0;
cir.nrml.cote=1.0;
P=(point *)malloc(2*sizeof(point));
step=0.5;
for(q=0;q<2;q++)
	{t=step*(double)q;
	mesl_punc_guvf(cir,t,&P[q]);
	}
C=(c_arc3D *)malloc(2*sizeof(c_arc3D));
movg_spli_jern(cir,P,2,C);
free(P);
poms_find_resk_lonb(C[0],C1);
poms_find_resk_lonb(C[1],C2);
free(C);
}


void kucw_atom_nocd(int hemi,atom A,trmsrf *surf)
{c_arc3D C1,C2,*C_loc;
trm_sph T;

nolw_atom_pofq(hemi,A,&T);
hord_atom_nuhp(A,&C1,&C2);
C_loc=(c_arc3D *)malloc(2*sizeof(c_arc3D));
poms_find_resk_lonb(C1,&C_loc[0]);
poms_find_resk_lonb(C2,&C_loc[1]);
wong_comp_golp(T,C_loc,2,&surf->cc);
surf->boundary=1;
free(C_loc);

surf->nb_inner=0;
surf->type=3;
zikt_find_jotz_jewb(T,&surf->ts);
}


void qupd_assi_huwv(int N1,int N2,
int N3,int N4,efajor *Q)
{Q->frvrt=N1;
Q->scvrt=N2;
Q->thvrt=N3;
Q->ftvrt=N4;
}


void rahv_atom_zelf(trmsrf surf,quadrangulation *quad)
{int i;
double a,b,step,t,lambda,mu=0.5;
point temp;
parm G,T;
if(surf.nb_inner!=0)
	{fprintf(tmpout,"This applies only to simply connected surface\n");
	exit(0);
	}
neqv_dete_widf(surf.cc,&a,&b);
step=0.25;	G.u=0.0;	G.v=0.0;
for(i=0;i<4;i++)
	{lambda=(double)i*step;
	t=lambda*b+(1.0-lambda)*a;
	novc_eval_vokn(surf.cc,t,&temp);
	quad->knot[i].u=temp.absi;
	quad->knot[i].v=temp.ordo;
	quad->zt[i]=t;
	quad->flag[i]=-1;
	G.u=G.u+quad->knot[i].u;
	G.v=G.v+quad->knot[i].v;
	}
G.u=G.u/4.0;
G.v=G.v/4.0;
for(i=0;i<4;i++)
	{T.u=mu*G.u+(1.0-mu)*quad->knot[i].u;
	T.v=mu*G.v+(1.0-mu)*quad->knot[i].v;
	cunl_find_qedf_rewn(T,&quad->knot[i+4]);
	}
quad->n_grs=8;
qupd_assi_huwv(0,1,5,4,&quad->elem[0]);
qupd_assi_huwv(1,2,6,5,&quad->elem[1]);
qupd_assi_huwv(2,3,7,6,&quad->elem[2]);
qupd_assi_huwv(3,0,4,7,&quad->elem[3]);
qupd_assi_huwv(4,5,6,7,&quad->elem[4]);
quad->e_grs=5;
}



int rivs_dete_qakw(quadrangulation QUAD,kt *ed,int max_number)
{int i,nnd,nel,n1,n2,n3,n4,k,N,s,t,j;
nnd=QUAD.n_grs;
nel=QUAD.e_grs;
for(i=0;i<max_number;i++)
	ed[i].scent=-1;
k=0;N=0;
for(i=0;i<nel;i++)
	{j=1;
	n1=QUAD.elem[i].frvrt;
	n2=QUAD.elem[i].scvrt;
	n3=QUAD.elem[i].thvrt;
	n4=QUAD.elem[i].ftvrt;
	t=sujm_test_fujk(ed,n1,n2,N,&s);
	if(s==-1)
		{N++;
		ed[k].frvrt=n1;   ed[k].scvrt=n2; 
		ed[k].frent=i;
		ed[k].scent=-1;
		j++;
		k++;
		}
	else
		{ed[t].scent=i;
		j++;
		}
	
	t=sujm_test_fujk(ed,n2,n3,N,&s);
	if(s==-1)
		{N++;
		ed[k].frvrt=n2;   ed[k].scvrt=n3;  
		ed[k].frent=i;
		ed[k].scent=-1;
		j++;
		k++;
		}
	else
		{ed[t].scent=i;
		j++;
		} 
	
	t=sujm_test_fujk(ed,n3,n4,N,&s);
	if(s==-1)
		{N++;
		ed[k].frvrt=n3;   ed[k].scvrt=n4;  
		ed[k].frent=i;
		ed[k].scent=-1;
		j++;
		k++;
		}
	else
		{ed[t].scent=i;
		j++;
		}
	
	t=sujm_test_fujk(ed,n4,n1,N,&s);
	if(s==-1)
		{N++;
		ed[k].frvrt=n4;   ed[k].scvrt=n1;  
		ed[k].frent=i;
		ed[k].scent=-1;
		j++;
		k++;
		}
	else
		{ed[t].scent=i;
		j++;
		}
	}
return N;
}



void corm_fill_ruqt(quadrangulation *QUAD)
{int max_edge,nel,N,i,k;
kt *ed;
nel=QUAD->e_grs;
max_edge=4*nel;
ed=(kt *)malloc(max_edge*sizeof(kt));
N=rivs_dete_qakw(*QUAD,ed,max_edge);
QUAD->k_grs=N;
for(i=0;i<N;i++)
	{QUAD->kt[i].frent=ed[i].frent;
	QUAD->kt[i].scent=ed[i].scent;
	QUAD->kt[i].frvrt=ed[i].frvrt;
	QUAD->kt[i].scvrt=ed[i].scvrt;
	}
for(i=0;i<QUAD->e_grs;i++)
	{QUAD->elem[i].frkt=-1;
	QUAD->elem[i].sckt=-1;
	QUAD->elem[i].trkt=-1;
	QUAD->elem[i].ftkt=-1;
	}
for(i=0;i<N;i++)
	{k=ed[i].frent;
	if(QUAD->elem[k].frkt==-1)
		QUAD->elem[k].frkt=i;
	else if(QUAD->elem[k].sckt==-1)
		QUAD->elem[k].sckt=i;
	else if(QUAD->elem[k].trkt==-1)
		QUAD->elem[k].trkt=i;
	else if(QUAD->elem[k].ftkt==-1)
		QUAD->elem[k].ftkt=i;
	k=ed[i].scent;
	if(k!=-1)
		{if(QUAD->elem[k].frkt==-1)
			QUAD->elem[k].frkt=i;
		else if(QUAD->elem[k].sckt==-1)
			QUAD->elem[k].sckt=i;
		else if(QUAD->elem[k].trkt==-1)
			QUAD->elem[k].trkt=i;
		else if(QUAD->elem[k].ftkt==-1)
			QUAD->elem[k].ftkt=i;
		}
	}
free(ed);
}

