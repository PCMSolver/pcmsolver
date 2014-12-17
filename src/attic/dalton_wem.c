#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include "constants.h"
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "sparse2.h"
#include "basis.h"
#include "dwt.h"
#include "WEM.h"
#include "WEMRHS.h"
#include "WEMPCG.h"
#include "WEMPGMRES.h"
#include "compression.h"
#include "init_points.h"
#include "interpolate.h"
#include "read_points.h"
#include "postproc.h"
#include "topology.h"
#include "energy.h"
#include "cubature.h"
#include "data.h"
#include "kern.h"
#include "gauss_square.h"
#include "dalton_wem.h"


#if !defined pi
  #define pi 3.1415926535897932385
#endif

sparse2         S;	      	/* komprimierte Steifigkeitsmatrix      */
vector3		*P;		/* Punkteliste der Einskalenbasis       */
unsigned int	**F;		/* Elementliste der Einskalenbasis      */
vector3		****T;		/* Oberflaecheninterpolation            */
vector3		***U;		/* Knotenpunkte                         */
unsigned int	nf;           	/* Laenge von F  	      	        */
unsigned int	p;		/* Anzahl der Patches                   */
unsigned int	M;	      	/* 2^M*2^M Elemente pro Patch           */
element		*E;		/* hierarchische Elementliste           */
wavelet		*W;		/* Liste der Wavelets                   */
double		*v;		/* approximierter Dichtevektor          */

/* 
 * Initializes the great WEM machine!
 *
 *   int (out) - number of potential points
 *
 * return 0 - if all ok
 *
 */
int dalton_wem_initialize_(int *num_points, double *da, double *db, double *ddp, 
			   double *apriori, double *aposteriori){

  char filename[80]="cavity.dat";

  *da=a;
  *db=b;
  *ddp=dp;

  /* Initialisierung */
  read_points(&U,&p,&M,filename);
  nf = p*(1<<M)*(1<<M);	        /* Anzahl der Patches */
  
  /* Topologie bestimmen */
  init_interpolate(&T,U,p,M);
  gennet(&P,&F,U,p,M);
  free_points(&U,p,M);

  /* erstelle Element-/Waveletliste */
  generate_elementlist(&E,P,F,p,M);
  generate_waveletlist(&W,E,p,M);
  set_quadrature_level(W,E,p,M);
  simplify_waveletlist(W,E,p,M);
  complete_elementlist(W,E,p,M);
  
  /* berechne a-priori Kompression */
  compression(&S,W,E,p,M,apriori);

  /* Steifigkeitsmatrix aufstellen */
  WEM(&S,W,P,E,T,p,M);

  /* berechne a-posteriori */
  postproc(&S,W,E,p,M,aposteriori);

  /* Number of potential points should be something like
   * [Number of parameter domains]*[2^(patch level)]**2*[quadrature level]**2
   */

  // Small chance for an error here..
  *num_points = nf * (general_quadrature_level+1)*(general_quadrature_level+1);

  v = NULL;

  if(general_quadrature_level!=0) {
    printf("Varo saatana! general_quadrature_level!=0!\n");
    return -1;
  }

  return 0; // perhaps one should check for unsufficient memory etc.
}

/*
 * Gets points for potential evaluation
 *
 *   data[number of potential points] (out) :
 *     [parameter domain][patch_level_y][patch_level_x][quadrature point]
 *   double * (out) - [x1,x2, ...]
 *   double * (out) - [y1,y2, ...]
 *   double * (out) - [z1,z2, ...]
 *     memory should be provided by dalton
 *
 */
void dalton_wem_get_potential_coordinates_(double *xp, double *yp, double *zp, double *carea){

  unsigned int i1, i2, i3, k, g, index;
  unsigned int n = 1 << M;
  vector2 s,t;
  vector3 res;
  double h = 1./n;
  cubature *Q;
  double totalArea=0.0;

  g = general_quadrature_level;
  init_Gauss_Square(&Q,g+1);

  index = 0;
  for (i1=0; i1<p; i1++){ /* sum over parameter domains */
    s.y = 0;
    for (i2=0; i2<n; i2++){
      s.x = 0;
      for (i3=0; i3<n; i3++){ 	/* sum over patches */
	
	for (k=0; k<Q[g].nop; k++){  /* sum over quadrature (g=2->nop=3*3) */
	  t = vector2_add(s,vector2_Smul(h,Q[g].xi[k]));
	  res = Chi(t,T[i1],M);
	  xp[index] = res.x;
	  yp[index] = res.y;
	  zp[index] = res.z;
	  totalArea += h*h*Q[g].w[k]*vector3_norm(n_Chi(t,T[i1],M));
	  index += 1;
	  }
	
	s.x += h;
      }
      s.y += h;
    }
  }
  
  *carea=totalArea;
}


void dalton_wem_get_areas_normals_(double *as, double *xn, double *yn, double *zn){

  unsigned int i1, i2, i3, k, g, index;
  unsigned int n = 1 << M;
  vector2 s,t;
  vector3 res;
  double h = 1./n;
  cubature *Q;
  //double totalArea=0.0;

  g = general_quadrature_level;
  init_Gauss_Square(&Q,g+1);

  index = 0;
  for (i1=0; i1<p; i1++){ /* sum over parameter domains */
    s.y = 0;
    for (i2=0; i2<n; i2++){
      s.x = 0;
      for (i3=0; i3<n; i3++){ 	/* sum over patches */
	
	for (k=0; k<Q[g].nop; k++){  /* sum over quadrature (g=2->nop=3*3) */
	  t = vector2_add(s,vector2_Smul(h,Q[g].xi[k]));
	  res = n_Chi(t,T[i1],M);
	  xn[index] = res.x;
	  yn[index] = res.y;
	  zn[index] = res.z;
	  as[index] = h*h*Q[g].w[k]*vector3_norm(n_Chi(t,T[i1],M));
	  index += 1;
	  }
	
	s.x += h;
      }
      s.y += h;
    }
  }
}


/*
 * Calculates surface polarization charges
 *
 *   int * (in)       - maximum number of charge points
 *   int * (out)      - number of charges
 *   double * (in)  - potential values at pre-specified points
 *   double * (out) - resulting charges at the same points
 *                    (allocate in DALTON, should be larger than #(potential))
 *   double * (out) - positions of charges [x1,x2, ...]
 *   double * (out) - positions of charges [y1,y2, ...]
 *   double * (out) - positions of charges [z1,z2, ...]
 *
 */
void dalton_wem_get_charges_(int *max_charges, int *num_charges, double *potential, double *charges){

  // This version uses only potential (not field)

  unsigned int i;

  // nf on oikea maara varauksille! 
  // potentiaalia voidaan evaluoida useammassakin paikkaa, mutta se on painotettu
  // quadraturen kertoimilla s.e. varaus on vakio (tai varaus*quadrature * potentiaal)

  *num_charges=nf*(general_quadrature_level+1)*(general_quadrature_level+1);
  if(*max_charges<*num_charges) {
    *num_charges=0;
    return;
  }

  if(v==NULL) v = (double*) calloc(*num_charges,sizeof(double)); /* initial guess? */
  memset(v,*num_charges*sizeof(double),0);
  memset(charges,*max_charges*sizeof(double),0);

  WEMRHS2D(charges,W,E,T,p,M,potential);
  WEMPCG(&S,charges,v,eps,p,M);

  tdwtKon(v,M,nf);

  double h=1.0/(1<<M);
  for(i=0;i<nf;i++) charges[i]=v[i]*h; // kertaa quadrature level... <- oikeasti kunnon quadrature kertoimilla
  // *num_charges=nf;
}

void dalton_wem_get_charges2_(int *max_charges, int *num_charges, double *fieldx,
			      double *fieldy,double *fieldz, double *charges){

  // Field version

  unsigned int i;

  // nf on oikea maara varauksille! 
  // potentiaalia voidaan evaluoida useammassakin paikkaa, mutta se on painotettu
  // quadraturen kertoimilla s.e. varaus on vakio (tai varaus*quadrature * potentiaal)

  *num_charges=nf*(general_quadrature_level+1)*(general_quadrature_level+1);
  if(*max_charges<*num_charges) {
    *num_charges=0;
    return;
  }

  if(v==NULL) v = (double*) calloc(*num_charges,sizeof(double)); /* initial guess? */

  WEMRHS1D(charges,W,E,T,p,M,fieldx,fieldy,fieldz);
  WEMPGMRES1(&S,charges,v,eps,p,M);

  tdwtKon(v,M,nf);

  double h=1.0/(1<<M);
  // kertaa quadrature level... <- oikeasti pitäis painottaa potentiaalia quadraturella (kai)
  for(i=0;i<nf;i++) charges[i]=v[i]*h;
  *num_charges=nf;
}


/*
 * Calculates the total energy
 *
 * double *potential;   Potential points!
 * double *energy;      Results
 *
 */
void dalton_wem_energy_(double *potential, double *energy){

  unsigned int	n = 1 << M;	/* n*n Patches pro Parametergebiet             */
  unsigned int	i1, i2, i3;	/* Laufindizes fuer Ansatzfunktion   	       */
  unsigned int	zi; 		/* Zeilenindex hieraus: zi = i1*(n*n)+i2*n+i3  */
  cubature      *Q;		/* Kubatur-Formeln                             */
  /*unsigned int	g = 2;	*/	/* Quadraturgrad                               */
  unsigned int	g = 0;		/* g=1, no difference wrt g=2 */
  unsigned int	k;	 	/* Laufindex fuer Quadratur		       */
  vector2	s;		/* Linker, unterer Eckpunkt des Patches zi     */
  vector2	t;		/* Auswertepunkte der Gaussquadratur	       */
  double	E;		/* Energie                                     */
  double	h = 1./n;	/* Schrittweite   			       */
  double	C = 0.0;
  unsigned int  index;

  if(v==NULL) exit(-102);

  /* Initialisierung */
  g = general_quadrature_level;
  init_Gauss_Square(&Q,g+1);	/* Kubatur-Formeln */
  
  /* Berechnung des Fehlers */
  index=0;
  E = zi = 0;
  for (i1=0; i1<p; i1++) /* sum over parameter domains */
    {  s.y = 0;
      for (i2=0; i2<n; i2++)
	{  s.x = 0;
	  for (i3=0; i3<n; i3++) 	/* zeilenweise Durchnumerierung der Patches zi = (i1,i2,i3) */
	    {  for (k=0; k<Q[g].nop; k++) /* sum over quadrature points (nop=3*3 for g=2) */
		{  t = vector2_add(s,vector2_Smul(h,Q[g].xi[k]));
		  /* E += Q[g].w[k]*v[zi]*f(Chi(t,T[i1],m)); */
		  E += Q[g].w[k]*v[zi]*potential[index];
		  index += 1;
		  C += Q[g].w[k]*v[zi];
		}
	      s.x += h;
	      zi++;
	    }
	  s.y += h;
	}
    }
  E = -0.5*h*E;	/* correct scaling */
  
  free_Gauss_Square(&Q,g+1);
  *energy=E;
}

void dalton_wem_finalize_(){
  if(v!=NULL) free(v);
  free_waveletlist(&W,p,M);
  free_elementlist(&E,p,M);
  free_interpolate(&T,p,M);
  free_patchlist(&F,nf);
  free_sparse2(&S);
  free(P);
}  

