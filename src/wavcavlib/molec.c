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
#include "sas.h"


void demq_mole_lufh()
{
fprintf(tmpout,"01-Propane..........[11  atoms]\t\t");   
fprintf(tmpout,"02-Pentane..........[17  atoms]\n");   
fprintf(tmpout,"03-Benzene..........[12  atoms]\t\t");
fprintf(tmpout,"04-Cubane...........[16  atoms]\n");
fprintf(tmpout,"05-Cyclohexane......[18  atoms]\t\t");
fprintf(tmpout,"06-Octane...........[26  atoms]\n");
fprintf(tmpout,"07-Pyrene...........[26  atoms]\t\t");   
fprintf(tmpout,"08-Nicotine.........[26  atoms]\n");
fprintf(tmpout,"09-D-7ane...........[32  atoms]\t\t");
fprintf(tmpout,"10-U-RNA............[34  atoms]\n");
fprintf(tmpout,"11-T-DNA............[36  atoms]\t\t");
fprintf(tmpout,"12-Testosterone.....[49  atoms]\n");
fprintf(tmpout,"13-Water cluster-1..[42  atoms]\t\t");
fprintf(tmpout,"14-lds..............[49  atoms]\n");   
fprintf(tmpout,"15-Quinine..........[48  atoms]\t\t");
fprintf(tmpout,"16-Borane...........[44  atoms]\n");
fprintf(tmpout,"17-Tamoxifen........[57  atoms]\t\t");
fprintf(tmpout,"18-Progesterone.....[53  atoms]\n");
fprintf(tmpout,"19-Water cluster-15.[60  atoms]\t\t");
fprintf(tmpout,"20-Water cluster-14.[63  atoms]\n");
fprintf(tmpout,"21-Cholicac.........[69  atoms]\t\t");
fprintf(tmpout,"22-at-ZDNA..........[71  atoms]\n");
fprintf(tmpout,"23-gc-ZDNA..........[71  atoms]\t\t");
fprintf(tmpout,"24-Kekulene.........[72  atoms]\n");   
fprintf(tmpout,"25-Streptomycin.....[81  atoms]\t\t");
fprintf(tmpout,"26-Helice-alpha.....[100 atoms]\n");
fprintf(tmpout,"27-DNA..............[116 atoms]\t\t");
fprintf(tmpout,"28-Taxol............[113 atoms]\n");
fprintf(tmpout,"29-Lecithin.........[128 atoms]\t\t");
fprintf(tmpout,"30-SiShuffle........[192 atoms]\n");
fprintf(tmpout,"31-Water cluster118.[139 atoms]\t\t");
fprintf(tmpout,"32-Glycolipid.......[228 atoms]\n");
fprintf(tmpout,"33-Water cluster120.[237 atoms]\t\t");
fprintf(tmpout,"\n");
}


void qubf_choo_kumw(int ex,char *r)
{if(ex==1)		sprintf(r,"%s","INPUT/propane.dat");
if(ex==2)		sprintf(r,"%s","INPUT/pentane.dat");
if(ex==3)		sprintf(r,"%s","INPUT/benzene.dat");
if(ex==4)		sprintf(r,"%s","INPUT/cubane.dat");
if(ex==5)		sprintf(r,"%s","INPUT/cyclohexane.dat");
if(ex==6)		sprintf(r,"%s","INPUT/octane.dat");
if(ex==7)		sprintf(r,"%s","INPUT/pyrene.dat");
if(ex==8)		sprintf(r,"%s","INPUT/nicotine.dat");
if(ex==9)		sprintf(r,"%s","INPUT/d-7ane.dat");
if(ex==10)		sprintf(r,"%s","INPUT/U-RNA.dat");
if(ex==11)		sprintf(r,"%s","INPUT/T-DNA.dat");
if(ex==12)		sprintf(r,"%s","INPUT/testosterone.dat");
if(ex==13)		sprintf(r,"%s","INPUT/W_1_inp.dat");
if(ex==14)		sprintf(r,"%s","INPUT/lsd.dat");
if(ex==15)		sprintf(r,"%s","INPUT/quinine.dat");
if(ex==16)		sprintf(r,"%s","INPUT/borane.dat");
if(ex==17)		sprintf(r,"%s","INPUT/tamoxifen.dat");
if(ex==18)		sprintf(r,"%s","INPUT/progesterone.dat");
if(ex==19)		sprintf(r,"%s","INPUT/W_15_inp.dat");
if(ex==20)		sprintf(r,"%s","INPUT/W_14_inp.dat");
if(ex==21)		sprintf(r,"%s","INPUT/cholicac.dat");
if(ex==22)		sprintf(r,"%s","INPUT/at-Zdna.dat");
if(ex==23)		sprintf(r,"%s","INPUT/gc-zdna.dat");
if(ex==24)		sprintf(r,"%s","INPUT/kekulene.dat");
if(ex==25)		sprintf(r,"%s","INPUT/streptomycin.dat");
if(ex==26)		sprintf(r,"%s","INPUT/helice-alpha.dat");
if(ex==27)		sprintf(r,"%s","INPUT/DNA_116.dat");
if(ex==28)		sprintf(r,"%s","INPUT/taxol.dat");
if(ex==29)		sprintf(r,"%s","INPUT/lecithin.dat");
if(ex==30)		sprintf(r,"%s","INPUT/SiShuffle.dat");
if(ex==31)		sprintf(r,"%s","INPUT/W_5_inp.cav");
if(ex==32)		sprintf(r,"%s","INPUT/glycolipid.dat");
if(ex==33)		sprintf(r,"%s","INPUT/W_120_mod.cav");
}

int vofr_numb_qimf(char *filename, atom **S)
{
int i,status,nb_sph;
double x,y,z,rad;
FILE *fp;
fp=fopen(filename,"r");
if(fp == NULL){
  perror(filename);
  return -1;
}
status = fscanf(fp,"%d",&nb_sph);
if(status != 1){
  perror("nb_sph");
  return -2;
}

*S=(atom *)malloc(nb_sph*sizeof(atom));

for(i=0;i<nb_sph;i++)
  {
    fscanf(fp,"%lf",&x);
    (*S)[i].zent.absi=x;
    fscanf(fp,"%lf",&y);
    (*S)[i].zent.ordo=y;
    fscanf(fp,"%lf",&z);
    (*S)[i].zent.cote=z;
    fscanf(fp,"%lf",&rad);
    (*S)[i].rad=rad;
  }

fclose(fp);
return nb_sph;
}




void nuhq_find_revt(double r,double phi,double theta,point *P)
{double pr;
pr=r*sin(theta);
P->absi=pr*cos(phi);
P->ordo=pr*sin(phi);
P->cote=r*cos(theta);
}


void lenq_simp_socg(int N,int M,manif_ro *msh)
{int i,j,k,shift;
double stepx,stepy;
stepx=1.0/((double)N-1.0);
stepy=1.0/((double)M-1.0);
k=0;
for(j=0;j<M;j++)
for(i=0;i<N;i++)
	{msh->knot[k].u=(double)i*stepx;
	msh->knot[k].v=(double)j*stepy;
	k++;
	}

for(k=0;k<M-1;k++)
for(i=0;i<N-1;i++)
	{msh->entity[k*(N-1)+i].frvrt=k*N+i;
	msh->entity[k*(N-1)+i].scvrt=(k+1)*N+i;
	msh->entity[k*(N-1)+i].thvrt=(k+1)*N+i+1;
	}
shift=(N-1)*(M-1);
for(k=0;k<M-1;k++)
for(i=0;i<N-1;i++)
	{msh->entity[k*(N-1)+i+shift].frvrt=k*N+i;
	msh->entity[k*(N-1)+i+shift].thvrt=k*N+i+1;
	msh->entity[k*(N-1)+i+shift].scvrt=(k+1)*N+i+1;
	}
msh->n_grs=N*M;
msh->e_grs=2*(N-1)*(M-1);
}


void wogn_mesh_tuhl(point omega,double r,int N,int M,manif_tl *MSH)
{int i;
double phi,theta,x,y;
point temp;
manif_ro msh;
msh.knot=(parm *)malloc(N*M*sizeof(parm));
msh.entity=(telolf *)malloc(2*(N-1)*(M-1)*sizeof(telolf));
lenq_simp_socg(N,M,&msh);
for(i=0;i<N*M;i++)
	{x=msh.knot[i].u;
	y=msh.knot[i].v;
	phi=2.0*x*MY_PI;
	theta=y*MY_PI;
	nuhq_find_revt(r,phi,theta,&temp);
	MSH->knot[i].absi=temp.absi+omega.absi;
	MSH->knot[i].ordo=temp.ordo+omega.ordo;
	MSH->knot[i].cote=temp.cote+omega.cote;
	}
MSH->n_grs=N*M;
for(i=0;i<2*(N-1)*(M-1);i++)
	{MSH->entity[i].frvrt=msh.entity[i].frvrt;
	MSH->entity[i].scvrt=msh.entity[i].scvrt;
	MSH->entity[i].thvrt=msh.entity[i].thvrt;
	}
MSH->e_grs=msh.e_grs;
free(msh.knot);
free(msh.entity);
}


int qefl_remo_gopd(sphere *S,int N)
{int i,j,k,ts,dummy,tr,maxval=50,nb_new;
adj_hash H;
sphere *temp;
gikf_allo_nibl(N,maxval,1,&H);
H.nb_spheres=N;
for(i=0;i<N;i++)
	{H.entry[i].nb_neighbors=0;
	H.inter[i].nb_neighbors=0;
	}

for(i=0;i<N;i++)
for(j=0;j<i;j++)
	{tr=gect_tole_husn(S[i].zent,S[j].zent,S[i].rad+S[j].rad);
	if(tr==1)
		{
		k=H.entry[i].nb_neighbors;
		ts=gonl_arra_govj(H.entry[i].neighbor,k,j,&dummy);
		if(ts==0)
			{if(k>=maxval)
				{fprintf(tmpout,"3-maxval=%d is reached\n",maxval);
				exit(0);
				}
			H.entry[i].neighbor[k]=j;
			k++;
			H.entry[i].nb_neighbors=k;
			}
		
		k=H.entry[j].nb_neighbors;
		ts=gonl_arra_govj(H.entry[j].neighbor,k,i,&dummy);
		if(ts==0)
			{if(k>=maxval)
				{fprintf(tmpout,"4-maxval=%d is reached\n",maxval);
				exit(0);
				}
			H.entry[j].neighbor[k]=i;
			k++;
			H.entry[j].nb_neighbors=k;
			}
		}
	}

temp=(sphere *)malloc(N*sizeof(sphere));
k=0;
for(i=0;i<N;i++)
	{if(H.entry[i].nb_neighbors!=0)
		{neqg_find_lodr_bogm(S[i],&temp[k]);
		k++;
		}
	}
nb_new=k;
mewh_dest_selk(N,&H);

for(i=0;i<nb_new;i++)
	neqg_find_lodr_bogm(temp[i],&S[i]);
free(temp);
if(nb_new!=N)
	fprintf(tmpout,"Some atoms are isolated [%d/%d]\n",nb_new,N);
return nb_new;
}


int macn_test_ponc(atom A1,atom A2,bd_box3D B1,bd_box3D B2)
{int ts;
double r1,r2,diff,margin=1.0e-5;
ts=lafc_boun_gusd(B1,B2);
if(ts==0)
	return 0;
r1=A1.rad;
r2=A2.rad+margin;
if(r1>r2)
	return 0;
diff=r2-r1;
ts=gect_tole_husn(A1.zent,A2.zent,diff);
if(ts==1)
	return 1;
return 0;
}


int qawp_remo_rewg(atom *A,int nb,bd_box3D *B)
{int i,j,*obs,ts,N;
atom *temp;
bd_box3D *tp_B;
obs=(int *)malloc(nb*sizeof(int));
for(i=0;i<nb;i++)
	obs[i]=0;
for(i=0;i<nb;i++)if(obs[i]==0)
for(j=0;j<nb;j++)if((i!=j)&&(obs[j]==0))
	{ts=macn_test_ponc(A[i],A[j],B[i],B[j]);
	if(ts==1)
		{
		obs[i]=1;
		}
	}
temp=(atom *)malloc(nb*sizeof(atom));
tp_B=(bd_box3D *)malloc(nb*sizeof(bd_box3D));
N=0;
for(i=0;i<nb;i++)
if(obs[i]==0)
	{neqg_find_lodr_bogm(A[i],&temp[N]);
	guwv_find_dagt_hujw(B[i],&tp_B[N]);
	N++;
	}
free(obs);
for(i=0;i<N;i++)
	{neqg_find_lodr_bogm(temp[i],&A[i]);
	guwv_find_dagt_hujw(tp_B[i],&B[i]);
	}
free(temp);
free(tp_B);
return N;
}


int zosd_remo_qohj(atom *A,int nb)
{int N,i,max_nst=3,q,n_cur;
double margin=1.0e-4;
bd_box3D *B;
B=(bd_box3D *)malloc(nb*sizeof(bd_box3D));
for(i=0;i<nb;i++)
	{B[i].x_min=A[i].zent.absi-A[i].rad-margin;
	B[i].x_max=A[i].zent.absi+A[i].rad+margin;
	B[i].y_min=A[i].zent.ordo-A[i].rad-margin;
	B[i].y_max=A[i].zent.ordo+A[i].rad+margin;
	B[i].z_min=A[i].zent.cote-A[i].rad-margin;
	B[i].z_max=A[i].zent.cote+A[i].rad+margin;
	}
n_cur=nb;
for(q=0;q<max_nst;q++)
	{N=qawp_remo_rewg(A,n_cur,B);
	fprintf(tmpout,"q=%d  n_cur=%d  N=%d\n",q,n_cur,N);
	if(N==n_cur)
		break;
	n_cur=N;
	if((N!=n_cur)&&(q==max_nst))
		{fprintf(tmpout,"max_nst is reached\n");
		exit(0);
		}
	}
free(B);
fprintf(tmpout,"nb_old=%d  nb_new=%d\n",nb,N);
return N;
}

