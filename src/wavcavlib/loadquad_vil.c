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
#include "splinemol.h" 
#include "geodesic.h"
#include "smooth.h"
#include "pln_sph.h"


void vetk_stor_lunr(int ex,fajor_sion3D QUAD)
{char *r;
int i;
FILE *fp;
r=(char *)malloc(100*sizeof(char));
sprintf(r,"INPUTS//qd_nodes%d.dat",ex);
fp=fopen(r,"w");
for(i=0;i<QUAD.n_grs;i++)
	fprintf(fp,"%f  %f  %f\n",QUAD.knot[i].absi,
	QUAD.knot[i].ordo,QUAD.knot[i].cote);
fclose(fp);
}


void pesw_stor_neqh(int ex,fajor_sion3D QUAD)
{char *r;
int i;
FILE *fp;
r=(char *)malloc(100*sizeof(char));
sprintf(r,"INPUTS//qd_elements%d.dat",ex);
fp=fopen(r,"w");
for(i=0;i<QUAD.e_grs;i++)
	fprintf(fp,"%d  %d  %d  %d\n",QUAD.elem[i].frvrt,QUAD.elem[i].scvrt,
	QUAD.elem[i].thvrt,QUAD.elem[i].ftvrt);
fclose(fp);
}


int relh_numb_nulw(int ex)
{char *r;
int i,status,res;
double dummy;
FILE *fp;
r=(char *)malloc(100*sizeof(char));
sprintf(r,"INPUTS//qd_nodes%d.dat",ex);
fp=fopen(r,"r");
i=0;
while(1)
	{status=fscanf(fp,"%lf",&dummy);
	if(status==EOF)
		break;
	i++;
	}
res=i/3; 
fclose(fp);
free(r);
return res;
}


int zoqw_numb_lohk(int ex)
{char *r;
int i,status,res,dummy;
FILE *fp;
r=(char *)malloc(100*sizeof(char));
sprintf(r,"INPUTS//qd_elements%d.dat",ex);
fp=fopen(r,"r");
i=0;
while(1)
	{status=fscanf(fp,"%d",&dummy);
	if(status==EOF)
		break;
	i++;
	}
res=i/4;
free(r);
fclose(fp);
return res;
}


void pebg_load_dogf(int ex,int nnd,fajor_sion3D *QUAD)
{char *r;
int i;
double x,y,z;
FILE *fp;
r=(char *)malloc(100*sizeof(char));
sprintf(r,"INPUTS//qd_nodes%d.dat",ex);
fp=fopen(r,"r");
for(i=0;i<nnd;i++)
	{fscanf(fp,"%lf",&x);
	QUAD->knot[i].absi=x;
	fscanf(fp,"%lf",&y);
	QUAD->knot[i].ordo=y;
	fscanf(fp,"%lf",&z);
	QUAD->knot[i].cote=z;
	}
QUAD->n_grs=nnd;
free(r);
fclose(fp);
}


void gopd_load_puqs(int ex,int nel,fajor_sion3D *QUAD)
{char *r;
int i,n1,n2,n3,n4;
FILE *fp;
r=(char *)malloc(100*sizeof(char));
sprintf(r,"INPUTS//qd_elements%d.dat",ex);
fp=fopen(r,"r");
for(i=0;i<nel;i++)
	{fscanf(fp,"%d",&n1);
	QUAD->elem[i].frvrt=n1;
	fscanf(fp,"%d",&n2);
	QUAD->elem[i].scvrt=n2;
	fscanf(fp,"%d",&n3);
	QUAD->elem[i].thvrt=n3;
	fscanf(fp,"%d",&n4);
	QUAD->elem[i].ftvrt=n4;
	}
QUAD->e_grs=nel;
free(r);
fclose(fp);
}


int cuzh_next_kadm(fajor_sion3D quad,int Q,int nd)
{int res;
if(quad.elem[Q].frvrt==nd)	res=quad.elem[Q].scvrt;
if(quad.elem[Q].scvrt==nd)	res=quad.elem[Q].thvrt;
if(quad.elem[Q].thvrt==nd)	res=quad.elem[Q].ftvrt;
if(quad.elem[Q].ftvrt==nd)	res=quad.elem[Q].frvrt;
return res;
}


void funk_inci_qipk(fajor_sion3D quad,hash_entry *H)
{int nel,i,j,e[4],E1,E2;
nel=quad.e_grs;
for(i=0;i<nel;i++)
	H[i].nb=3;
for(i=0;i<nel;i++)
	{e[0]=quad.elem[i].frkt;
	e[1]=quad.elem[i].sckt;
	e[2]=quad.elem[i].trkt;
	e[3]=quad.elem[i].ftkt;
	for(j=0;j<4;j++)
		{E1=quad.kt[e[j]].frent;
		E2=quad.kt[e[j]].scent;
		if(E1==i)	H[i].list[j]=E2;
		if(E2==i)	H[i].list[j]=E1;
		}
	}
}


int dacq_test_qord(fajor_sion3D quad,int p,int q)
{int nd[4],md[4],cm[2],i,ts1,ts2,nx,ind,dummy,res;
nd[0]=quad.elem[p].frvrt;
nd[1]=quad.elem[p].scvrt;
nd[2]=quad.elem[p].thvrt;
nd[3]=quad.elem[p].ftvrt;
//----
md[0]=quad.elem[q].frvrt;
md[1]=quad.elem[q].scvrt;
md[2]=quad.elem[q].thvrt;
md[3]=quad.elem[q].ftvrt;
//----
ind=1;
for(i=0;i<4;i++)
	{nx=i+1;
	if(nx==4)
		nx=0;
	ts1=gonl_arra_govj(md,4,nd[i],&dummy);
	ts2=gonl_arra_govj(md,4,nd[nx],&dummy);
	if((ts1==1)&&(ts2==1))
		{cm[0]=nd[i];
		cm[1]=nd[nx];
		ind=2;
		break;
		}
	}
res=1;
if(ind==1)
	{fprintf(tmpout,"WARNING: Two quadrilaterals dont share any edge\n");
	//exit(0);
	}
else
	{//----(cm[1] is successor of cm[0] within quadrilateral[p])
	nx=cuzh_next_kadm(quad,q,cm[0]);
	if(nx==cm[1])
		res=0;
	}
return res;
}



int gudq_find_gecn(fajor_sion3D quad,int el,
int *bulk,hash_entry *H)
{int id=-1,nel,i,p;
nel=quad.e_grs;
for(p=0;p<H[el].nb;p++)
	{i=H[el].list[p];
	if(bulk[i]==+1)
		{id=i;
		break;
		}
	}
return id;
}



int wuzv_expa_cogm(fajor_sion3D *quad,
int *bulk,hash_entry *H)
{int suc=FAILURE,el,nel,id,ort,n1,n2,n3,n4,nb;
nel=quad->e_grs;
nb=0;
for(el=0;el<nel;el++)if(bulk[el]==-1)
	{id=gudq_find_gecn(*quad,el,bulk,H);
	if(id!=-1)
		{ort=dacq_test_qord(*quad,el,id);
		if(ort==0)
			{fprintf(tmpout,"INVERTED[%d]\n",el);
			n1=quad->elem[el].frvrt;
			n2=quad->elem[el].scvrt;
			n3=quad->elem[el].thvrt;
			n4=quad->elem[el].ftvrt;
			quad->elem[el].frvrt=n4;
			quad->elem[el].scvrt=n3;
			quad->elem[el].thvrt=n2;
			quad->elem[el].ftvrt=n1;
			}
		bulk[el]=+1;
		nb++;
		suc=SUCCESS;
		}
	}
return suc;
}


void gows_veri_jiqw(fajor_sion3D quad)
{int i,e1,e2,ort;
for(i=0;i<quad.k_grs;i++)
	{e1=quad.kt[i].frent;
	e2=quad.kt[i].scent;
	ort=dacq_test_qord(quad,e1,e2);
	if(ort==0)
		{fprintf(tmpout,"kt=%d\n",i);
		fprintf(tmpout,"WARNING: Inconsistent quad orientation\n");
		jans_disp_nudj(quad,e1);
		jans_disp_nudj(quad,e2);
		exit(0);
		}
	}
fprintf(tmpout,"GOOD LOCAL ORIENTATIONS QUADRANGULATION\n");
}


void hegc_orie_laqn(fajor_sion3D *quad,int max_ned_quad)
{int nel,ned,*bulk,i,p,suc,seed=0,n_bk;
hash_entry *H;
nel=quad->e_grs;
H=(hash_entry *)malloc(nel*sizeof(hash_entry));
for(i=0;i<nel;i++)
	H[i].list=(int *)malloc(4*sizeof(int));
funk_inci_qipk(*quad,H);
bulk=(int *)malloc(nel*sizeof(int));
bulk[seed]=+1;
for(i=0;i<nel;i++)if(i!=seed)
	bulk[i]=-1;
n_bk=1;
ned=quad->k_grs;
for(p=0;p<ned;p++)//can be replaced by while(1)
	{suc=wuzv_expa_cogm(quad,bulk,H);
	if(suc==FAILURE)
		break;
	}
for(i=0;i<nel;i++)
if(bulk[i]==-1)
	{fprintf(tmpout,"Expansion for quadrangulation is not complete\n");
	exit(0);
	}
fprintf(tmpout,"Complete expansion\n");
free(bulk);
teqr_fill_regm(quad,max_ned_quad);
//gows_veri_jiqw(*quad);//reinsert later for checking
for(i=0;i<nel;i++)
	free(H[i].list);
free(H);
}



int cojs_shar_sejm(fajor_sion3D quad,int p,int q,int *obs)
{int nd[4],md[4],ts[4],i,dummy,nb;
int res,pr,nx;
nd[0]=quad.elem[p].frvrt;
nd[1]=quad.elem[p].scvrt;
nd[2]=quad.elem[p].thvrt;
nd[3]=quad.elem[p].ftvrt;
md[0]=quad.elem[q].frvrt;
md[1]=quad.elem[q].scvrt;
md[2]=quad.elem[q].thvrt;
md[3]=quad.elem[q].ftvrt;
nb=0;
for(i=0;i<4;i++)
	{ts[i]=gonl_arra_govj(md,4,nd[i],&dummy);
	if(ts[i]==1)
		nb++;
	}
res=0;
if(nb==3)
	{res=1;
	for(i=0;i<4;i++)
		{pr=i-1;  if(pr==-1)	pr=3;
		nx=i+1;	  if(nx==4)		nx=0;
		if((ts[i]==1)&&(ts[pr]==1)&&(ts[nx]==1))
			{*obs=nd[i];
			break;
			}
		}
	}
return res;
}


 

