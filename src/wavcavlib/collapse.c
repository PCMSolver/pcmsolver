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
#include "splinemol.h"
#include "pln_sph.h"
#include "triang.h"
#include "coarsequad.h"


void roph_coll_fact(manif_tl *msh,int n1,int n2,point P)
{msh->knot[n1].absi=P.absi;
msh->knot[n1].ordo=P.ordo;
msh->knot[n1].cote=P.cote;

msh->knot[n2].absi=P.absi;
msh->knot[n2].ordo=P.ordo;
msh->knot[n2].cote=P.cote;
}


void tedn_find_bals_qucl(teboka_topo x,teboka_topo *y)
{int i,val;
val=x.val;
for(i=0;i<val;i++)
	y->inc[i]=x.inc[i];
y->val=x.val;
}



void gohl_disc_gunp(manif_tl *msh,int m)
{int nnd,nel,ned,*idx,i,k,n1,n2,n3;
point *temp;
teboka_topo *tp;
nnd =msh->n_grs;
idx =(int *)malloc(nnd*sizeof(int));
temp=(point *)malloc(nnd*sizeof(point));
tp  =(teboka_topo *)malloc(nnd*sizeof(teboka_topo));
k=0;
for(i=0;i<nnd;i++)
if(i!=m)
	{temp[k].absi=msh->knot[i].absi;
	temp[k].ordo=msh->knot[i].ordo;
	temp[k].cote=msh->knot[i].cote;
	tedn_find_bals_qucl(msh->increl[i],&tp[k]);
	idx[i]=k;
	k++;
	}
idx[m]=-1;

for(i=0;i<nnd-1;i++)
	{msh->knot[i].absi=temp[i].absi;
	msh->knot[i].ordo=temp[i].ordo;
	msh->knot[i].cote=temp[i].cote;
	tedn_find_bals_qucl(tp[i],&msh->increl[i]);
	}
msh->n_grs=nnd-1;

nel=msh->e_grs;
for(i=0;i<nel;i++)
	{n1=msh->entity[i].frvrt;	msh->entity[i].frvrt=idx[n1];
	n2=msh->entity[i].scvrt;		msh->entity[i].scvrt=idx[n2];
	n3=msh->entity[i].thvrt;		msh->entity[i].thvrt=idx[n3];
	}

ned=msh->k_grs;
for(i=0;i<ned;i++)
	{n1=msh->kt[i].frvrt;		msh->kt[i].frvrt=idx[n1];
	n2=msh->kt[i].scvrt;		msh->kt[i].scvrt=idx[n2];
	}
free(idx);
free(temp);
free(tp);
}



void fuwr_repl_wejf(manif_tl *msh,int n1,int n2)
{int i,k,nel,ned,na,nb,nc,*temp1,*temp2;
int val1,val2,*temp,dummy,ts,m1,m2,m3;
nel=msh->e_grs;
ned=msh->k_grs;
for(i=0;i<nel;i++)
	{na=msh->entity[i].frvrt;	if(na==n2)	msh->entity[i].frvrt=n1;
	nb=msh->entity[i].scvrt;		if(nb==n2)	msh->entity[i].scvrt=n1;
	nc=msh->entity[i].thvrt;		if(nc==n2)	msh->entity[i].thvrt=n1;
	m1=msh->entity[i].frvrt;
	m2=msh->entity[i].scvrt;
	m3=msh->entity[i].thvrt;
	if((m1==n1)||(m2==n1)||(m3==n1))
		msh->entity[i].ms_wert=valp_area_qelk(msh->knot[m1],msh->knot[m2],msh->knot[m3]);
	}
for(i=0;i<ned;i++)
	{na=msh->kt[i].frvrt;			if(na==n2)	msh->kt[i].frvrt=n1;
	nb=msh->kt[i].scvrt;			if(nb==n2)	msh->kt[i].scvrt=n1;
	}

val1=msh->increl[n1].val;
val2=msh->increl[n2].val;
temp1=(int *)malloc(val1*sizeof(int));
temp2=(int *)malloc(val2*sizeof(int));
for(i=0;i<val1;i++)
	temp1[i]=msh->increl[n1].inc[i];
for(i=0;i<val2;i++)
	temp2[i]=msh->increl[n2].inc[i];
temp=(int *)malloc((val1+val2)*sizeof(int));
k=0;
for(i=0;i<val1;i++)
	{ts=gonl_arra_govj(temp,k,temp1[i],&dummy);
	if((ts==0)&&(temp1[i]!=-1))
		{temp[k]=temp1[i];
		k++;
		}
	}
for(i=0;i<val2;i++)
	{ts=gonl_arra_govj(temp,k,temp2[i],&dummy);
	if((ts==0)&&(temp2[i]!=-1))
		{temp[k]=temp2[i];
		k++;
		}
	}
free(temp1);
free(temp2);
msh->increl[n1].val=k;
for(i=0;i<k;i++)
	msh->increl[n1].inc[i]=temp[i];
free(temp);
}



void zudp_disc_duwc(manif_tl *msh,int m,int *map)
{int i,k,nnd;
nnd=msh->n_grs;
gohl_disc_gunp(msh,m);
map[m]=-1;
k=0;
for(i=0;i<nnd;i++)
if(i!=m)
	{map[i]=k;
	k++;
	}
}



void vegs_disc_lepq(manif_tl *msh,int *el,int M)
{int i,k,nel,ts,dummy,*map,e1,e2;
telolf *temp;
nel=msh->e_grs;
temp=(telolf *)malloc((nel-M)*sizeof(telolf));
map=(int *)malloc(nel*sizeof(int));
k=0;
for(i=0;i<nel;i++)
	{ts=gonl_arra_govj(el,M,i,&dummy);
	if(ts==0)
		{temp[k].frvrt=msh->entity[i].frvrt;
		temp[k].scvrt=msh->entity[i].scvrt;
		temp[k].thvrt=msh->entity[i].thvrt;
		temp[k].frkt=msh->entity[i].frkt;
		temp[k].sckt=msh->entity[i].sckt;
		temp[k].trkt=msh->entity[i].trkt;
		temp[k].ms_wert=msh->entity[i].ms_wert;
		map[i]=k;
		k++;
		}
	else
		map[i]=-1;
	}

msh->e_grs=nel-M;
for(i=0;i<nel-M;i++)
	{msh->entity[i].frvrt=temp[i].frvrt;
	msh->entity[i].scvrt=temp[i].scvrt;
	msh->entity[i].thvrt=temp[i].thvrt;
	msh->entity[i].frkt=temp[i].frkt;
	msh->entity[i].sckt=temp[i].sckt;
	msh->entity[i].trkt=temp[i].trkt;
	msh->entity[i].ms_wert=temp[i].ms_wert;
	}
free(temp);

for(i=0;i<msh->k_grs;i++)
	{e1=msh->kt[i].frent;
	e2=msh->kt[i].scent;
	msh->kt[i].frent=map[e1];
	msh->kt[i].scent=map[e2];
	}
free(map);
}



void topr_disc_nufr(manif_tl *msh,int e,int *map)
{int i,j,k,nnd,ned,nel,*idx,ed1,ed2,ed3,val,s,*tp,m;
kt *temp;
ned=msh->k_grs;
temp=(kt *)malloc(ned*sizeof(kt));
idx=(int *)malloc(ned*sizeof(int));
map[e]=-1;
idx[e]=-1;
k=0;
for(i=0;i<ned;i++)if(i!=e)
	{temp[k].frvrt=msh->kt[i].frvrt;
	temp[k].scvrt=msh->kt[i].scvrt;
	temp[k].frent=msh->kt[i].frent;
	temp[k].scent=msh->kt[i].scent;
	idx[i]=k;
	map[i]=k;
	k++;
	}

nel=msh->e_grs;
for(i=0;i<nel;i++)
	{ed1=msh->entity[i].frkt;	msh->entity[i].frkt=idx[ed1];
	ed2=msh->entity[i].sckt;	msh->entity[i].sckt=idx[ed2];
	ed3=msh->entity[i].trkt;	msh->entity[i].trkt=idx[ed3];
	}

nnd=msh->n_grs;
for(i=0;i<nnd;i++)
	{val=msh->increl[i].val;
	tp=(int *)malloc(val*sizeof(int));
	m=0;
	for(j=0;j<val;j++)
		{s=msh->increl[i].inc[j];
		if(idx[s]!=-1)
			{tp[m]=idx[s];
			m++;
			}
		}
	for(j=0;j<m;j++)
		msh->increl[i].inc[j]=tp[j];
	msh->increl[i].val=m;
	free(tp);
	}

for(i=0;i<ned-1;i++)
	{msh->kt[i].frvrt=temp[i].frvrt;
	msh->kt[i].scvrt=temp[i].scvrt;
	msh->kt[i].frent=temp[i].frent;
	msh->kt[i].scent=temp[i].scent;
	}
free(temp);
free(idx);
msh->k_grs=ned-1;
}



void nigf_disc_pojm(manif_tl *msh,int *e,int M,int *map)
{int i,j,k,m,nnd,ned,nel,*idx,ed1,ed2,ed3,val,s,ts,dummy,*tp;
kt *temp;
ned=msh->k_grs;
idx=(int *)malloc(ned*sizeof(int));
for(i=0;i<M;i++)
	{map[e[i]]=-1;
	idx[e[i]]=-1;
	}
temp=(kt *)malloc(ned*sizeof(kt));
k=0;
for(i=0;i<ned;i++)
	{ts=gonl_arra_govj(e,M,i,&dummy);
	if(ts==0)
		{temp[k].frvrt=msh->kt[i].frvrt;
		temp[k].scvrt=msh->kt[i].scvrt;
		temp[k].frent=msh->kt[i].frent;
		temp[k].scent=msh->kt[i].scent;
		idx[i]=k;
		map[i]=k;
		k++;
		}
	}

nel=msh->e_grs;
for(i=0;i<nel;i++)
	{ed1=msh->entity[i].frkt;	msh->entity[i].frkt=idx[ed1];
	ed2=msh->entity[i].sckt;	msh->entity[i].sckt=idx[ed2];
	ed3=msh->entity[i].trkt;	msh->entity[i].trkt=idx[ed3];
	}

nnd=msh->n_grs;
for(i=0;i<nnd;i++)
	{val=msh->increl[i].val;
	tp=(int *)malloc(val*sizeof(int));
	m=0;
	for(j=0;j<val;j++)
		{s=msh->increl[i].inc[j];
		if(idx[s]!=-1)
			{tp[m]=idx[s];
			if(idx[s]<0)
				{fprintf(tmpout,"inaccurate indices\n");
				exit(0);
				}
			m++;
			}
		}
	for(j=0;j<m;j++)
		msh->increl[i].inc[j]=tp[j];
	msh->increl[i].val=m;
	free(tp);
	}

for(i=0;i<ned-M;i++)
	{msh->kt[i].frvrt=temp[i].frvrt;
	msh->kt[i].scvrt=temp[i].scvrt;
	msh->kt[i].frent=temp[i].frent;
	msh->kt[i].scent=temp[i].scent;
	}
free(temp);
free(idx);
msh->k_grs=ned-M;
}



void dopr_fuse_tenv(manif_tl *msh,int e1,int e2,int el3,int el4,
int e3,int e4,int el5,int el6,int *map)
{int *sev,E1,E3;
if(msh->entity[el3].frkt==e2)	msh->entity[el3].frkt=e1;	
if(msh->entity[el3].sckt==e2)	msh->entity[el3].sckt=e1;		
if(msh->entity[el3].trkt==e2)	msh->entity[el3].trkt=e1;	

if(msh->entity[el4].frkt==e2)	msh->entity[el4].frkt=e1;	
if(msh->entity[el4].sckt==e2)	msh->entity[el4].sckt=e1;		
if(msh->entity[el4].trkt==e2)	msh->entity[el4].trkt=e1;	

if(msh->entity[el5].frkt==e4)	msh->entity[el5].frkt=e3;	
if(msh->entity[el5].sckt==e4)	msh->entity[el5].sckt=e3;		
if(msh->entity[el5].trkt==e4)	msh->entity[el5].trkt=e3;	

if(msh->entity[el6].frkt==e4)	msh->entity[el6].frkt=e3;	
if(msh->entity[el6].sckt==e4)	msh->entity[el6].sckt=e3;		
if(msh->entity[el6].trkt==e4)	msh->entity[el6].trkt=e3;	

sev=(int *)malloc(2*sizeof(int));
sev[0]=e2;
sev[1]=e4;
nigf_disc_pojm(msh,sev,2,map);
free(sev);
E1=map[e1];
E3=map[e3];
msh->kt[E1].frent=el3;	msh->kt[E1].scent=el4;
msh->kt[E3].frent=el5;	msh->kt[E3].scent=el6;
}


double juhz_area_hong(manif_tl msh,int s)
{int n1,n2,n3;
double surf;
n1=msh.entity[s].frvrt;
n2=msh.entity[s].scvrt;
n3=msh.entity[s].thvrt;
surf=valp_area_qelk(msh.knot[n1],msh.knot[n2],msh.knot[n3]);
return surf;
}



int sokm_loca_vazn(manif_tl msh,int s)
{int ts=1,*nd,ed[3],i,e,n1,n2,ned;
int tr1,tr2,dummy,nel;
double surf;
ned=msh.k_grs;
nd=(int *)malloc(3*sizeof(int));
nd[0]=msh.entity[s].frvrt;
nd[1]=msh.entity[s].scvrt;
nd[2]=msh.entity[s].thvrt;
ed[0]=msh.entity[s].frkt;
ed[1]=msh.entity[s].sckt;
ed[2]=msh.entity[s].trkt;
if((ed[0]<0)||(ed[1]<0)||(ed[2]<0))
	{fprintf(tmpout,"edges of el[%d]=[%d,%d,%d]\n",s,ed[0],ed[1],ed[2]);
	exit(0);
	}
if((ed[0]>=ned)||(ed[1]>=ned)||(ed[2]>=ned))
	{fprintf(tmpout,"edges of el[%d]=[%d,%d,%d]\n",s,ed[0],ed[1],ed[2]);
	fprintf(tmpout,"ned=%d\n",ned);
	exit(0);
	}

nel=msh.e_grs;
for(i=0;i<3;i++)
	{e=ed[i];
	n1=msh.kt[e].frvrt;
	n2=msh.kt[e].scvrt;
	tr1=gonl_arra_govj(nd,3,n1,&dummy);
	tr2=gonl_arra_govj(nd,3,n2,&dummy);
	if((tr1==0)||(tr2==0))
		{ts=0;
		break;
		}
	}
free(nd);
if(ts==0)
	{fprintf(tmpout,"element=%d with bad integrity1\n",s);
	exit(0);
	}

for(i=0;i<3;i++)
	{e=ed[i];
	if((msh.kt[e].frent!=s)&&(msh.kt[e].scent!=s))
		{ts=0;
		break;
		}
	}
if(ts==0)
	{fprintf(tmpout,"[nd1,nd2,nd3]=[%d,%d,%d]\n",
	msh.entity[s].frvrt,msh.entity[s].scvrt,msh.entity[s].thvrt);
	fprintf(tmpout,"[ed1,ed2,ed3]=[%d,%d,%d]\n",
	msh.entity[s].frkt,msh.entity[s].sckt,msh.entity[s].trkt);
	for(i=0;i<3;i++)
		fprintf(tmpout,"ed=%d    [n1,n2]=[%d,%d]   inc=[%d,%d]\n",
		ed[i],msh.kt[ed[i]].frvrt,msh.kt[ed[i]].scvrt,msh.kt[ed[i]].frent,msh.kt[ed[i]].scent);
	fprintf(tmpout,"element=%d with bad integrity2\n",s);
	exit(0);
	}

surf=juhz_area_hong(msh,s);
if(surf<1.0e-4)
	{ts=0;
	
	}
return ts;
}



int cugt_loca_pogt(manif_tl msh,int s)
{int ts=1,e[2],i,j,f,n[2],r,tr,val;
e[0]=msh.kt[s].frent;
e[1]=msh.kt[s].scent;

for(i=0;i<2;i++)
	{f=e[i];
	if((msh.entity[f].frkt!=s)&&(msh.entity[f].sckt!=s)&&(msh.entity[f].trkt!=s))
		{ts=0;
		break;
		}
	}
if(ts==0)
	{fprintf(tmpout,"kt=%d with bad integrity1\n",s);
	fprintf(tmpout,"element1=%d  incident edges=[%d,%d,%d]\n",
	e[0],msh.entity[e[0]].frkt,msh.entity[e[0]].sckt,msh.entity[e[0]].trkt);
	fprintf(tmpout,"element2=%d  incident edges=[%d,%d,%d]\n",
	e[1],msh.entity[e[1]].frkt,msh.entity[e[1]].sckt,msh.entity[e[1]].trkt);
	exit(0);
	}

n[0]=msh.kt[s].frvrt;
n[1]=msh.kt[s].scvrt;
for(i=0;i<2;i++)
	{r=n[i];
	val=msh.increl[r].val;
	tr=-1;
	for(j=0;j<val;j++)
		{if(msh.increl[r].inc[j]==s)
			{tr=+1;
			break;
			}
		}
	if(tr==-1)
		{fprintf(tmpout,"kt=%d with bad integrity2 with node=%d\n",s,n[i]);
		exit(0);	
		}
	}
return ts;
}


int darj_loca_wejf(manif_tl msh,int s)
{int ts=1,e[2],i,f;
e[0]=msh.kt[s].frent;
e[1]=msh.kt[s].scent;

for(i=0;i<2;i++)
	{f=e[i];
	if((msh.entity[f].frkt!=s)&&(msh.entity[f].sckt!=s)&&(msh.entity[f].trkt!=s))
		{ts=0;
		break;
		}
	}
if(ts==0)
	{fprintf(tmpout,"kt=%d with bad integrity1\n",s);
	fprintf(tmpout,"element1=%d  incident edges=[%d,%d,%d]\n",
	e[0],msh.entity[e[0]].frkt,msh.entity[e[0]].sckt,msh.entity[e[0]].trkt);
	fprintf(tmpout,"element2=%d  incident edges=[%d,%d,%d]\n",
	e[1],msh.entity[e[1]].frkt,msh.entity[e[1]].sckt,msh.entity[e[1]].trkt);
	exit(0);
	}
return ts;
}



