/************
 *  main.c  *
 ************/


/*=================================================*
 *  Hauptprogramm zum Wavelet-Galerkin-Verfahren.  *
 *=================================================*/
 
 
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "sparse2.h"
#include "sparse.h"
#include "basis.h"
#include "dwt.h"
#include "WEM.h"
#include "WEMRHS_pwl.h"
#include "WEMPCG_pwl.h"
#include "WEMPGMRES_pwl.h"
#include "compression.h"
#include "interpolate.h"
#include "read_points.h"
#include "postproc.h"
#include "topology.h"
#include "precond.h"
#include "energy.h"
#include "volume.h"
#include "kern.h"


#if !defined pi
  #define pi 3.1415926535897932385
#endif  


int main()
{
sparse2         S_i, S_e;	/* komprimierte Steifigkeitsmatrix    */
sparse		G;		/* Massenmatrix			      */
vector3		*P;		/* Punkteliste der Einskalenbasis     */
unsigned int	**F;		/* Elementliste der Einskalenbasis    */
vector3		****T;		/* Oberflaecheninterpolation          */
vector3		***U;		/* Knotenpunkte                       */
unsigned int	nf;           	/* Laenge von F  	      	      */
unsigned int	np;           	/* Laenge von P  	      	      */
unsigned int	p;		/* Anzahl der Patches                 */
unsigned int	M;	      	/* 2^M*2^M Elemente pro Patch         */
element		*E;		/* hierarchische Elementliste         */
wavelet		*W;		/* Liste der Wavelets                 */
double		*u, *v;		/* approximierter Dichtevektor        */
double		*rhs;		/* rechte Seite des Gleichungssystems */
double		res;		/* berechnete Austauschenergien       */ 
time_t		t1, t2, t3;	/* Sartzeit/Zwischenzeit/Endzeit      */
unsigned int	i, j;		/* Laufindizes                        */
double		det;		/* Determinante im anistropen Fall    */

/*========================================================*/
unsigned int	CASE = 1;	/* FLAG FOR CHOICE OF BIE */
/*========================================================*/

/* Ausgabe */
switch (CASE)
{  case 1:  printf("Pure Poisson Equation with 2nd kind integral equation:\n");
            printf("======================================================\n\n");
            break;
   case 2:  printf("Pure Poisson Equation with 1st kind integral equation:\n");
            printf("======================================================\n\n");
            break;
   case 3:  printf("Poisson-Boltzmann Equation:\n");
            printf("===========================\n\n");
            break;
   case 4:  printf("Anisotropic Dielectric:\n");
            printf("=======================\n\n");
            break;
   default: printf("ERROR: non-existent scheme\n"); return(0);
   }

/* Initialisierung */
read_points(&U,&p,&M);
nf = p*(1<<M)*(1<<M);	        /* Anzahl der Patches */
time(&t1);

/* Topologie bestimmen */
init_interpolate(&T,U,p,M);
np = gennet(&P,&F,U,p,M);
free_points(&U,p,M);
volume(F,T,p,M);

/* erstelle Element-/Waveletliste */
printf("Number of levels:                %d \n",M);
printf("Number of parameter domains:     %d \n",p);
printf("Number of ansatz functions:      %d \n",np);
printf("Computing the overhead:          ");
generate_elementlist(&E,P,F,p,M);
generate_waveletlist(&W,E,p,M,np);
set_quadrature_level(W,E,p,M,np);
simplify_waveletlist(W,E,p,M,np);
complete_elementlist(W,E,p,M,np);
time(&t2);
printf("Computation time:                %g secs.\n\n",difftime(t2,t1));

/* erstes Paar Steifigkeitsmatrizen aufstellen */
printf("Computing the 1st pair of system matrices: \n");
compression(&S_i,W,E,p,M,np);
if (CASE < 3) WEM_pwl(&S_i,W,P,E,T,p,M,SingleLayerInt,DoubleLayerInt,2*pi*(1+epsilon)/(1-epsilon));
else          WEM_pwl(&S_i,W,P,E,T,p,M,SingleLayerInt,DoubleLayerInt,2*pi);
postproc(&S_i,W,E,p,M);
time(&t3);      /* Zwischenzeit */
printf("Computation time:                %g secs.\n\n",difftime(t3,t2));

/* zweites Paar Steifigkeitsmatrizen aufstellen */
if ((CASE == 3) || (CASE == 4))
{  time(&t2);
   printf("Computing the 2nd pair of system matrices: \n");
   compression(&S_e,W,E,p,M,np);
   if (CASE == 3) 
   {  WEM_pwl(&S_e,W,P,E,T,p,M,SingleLayerExt,DoubleLayerExt,-2*pi);
      for (i=0; i<np; i++)	/* correct scaling */
      {  for (j=0; j<S_e.row_number[i]; j++) S_e.value1[i][j] /= epsilon;
         }
      }
   else
   {  det = sqrt( epsilon11*epsilon22*epsilon33 + epsilon12*epsilon23*epsilon31 + epsilon13*epsilon21*epsilon32 \
                - epsilon11*epsilon32*epsilon23 - epsilon21*epsilon12*epsilon33 - epsilon31*epsilon22*epsilon13 );   
      WEM_pwl(&S_e,W,P,E,T,p,M,SingleLayerAni,DoubleLayerAni,-2*pi/det);
      for (i=0; i<np; i++)	/* correct scaling */
      {  for (j=0; j<S_e.row_number[i]; j++) 
         {  S_e.value1[i][j] *= det;
            S_e.value2[i][j] *= det;
	    } 
         }
      }
   postproc(&S_e,W,E,p,M);
   time(&t3);      /* Zwischenzeit */
   printf("Computation time:                %g secs.\n\n",difftime(t3,t2));
   }

/* loese Gleichungssystem */
u = (double*) calloc(np,sizeof(double));
v = (double*) calloc(np,sizeof(double));
if (CASE == 1)
{  WEMRHS_pwl1(&rhs,W,E,T,p,M,np);
   i = WEMPGMRES_pwl1(&S_i,rhs,u,eps,W,F,p,M);
   printf("Solving the linear system:       %d iterations\n",i);
   }
else if (CASE == 2)
{  WEMRHS_pwl2(&rhs,W,E,T,p,M,np);		/* compute correct rhs: b-G*A2^(-1)*b */
   i = WEMPGMRES_pwl2(&S_i,rhs,v,eps,W,F,p,M);
   printf("Solving the 1st linear system:   %d iterations\n",i);
   init_sparse(&G,np,np,10);
   single_scale_gram(&G,F,p,M);
   tdwtLin(v,F,M,p,np);
   for (i=0; i<np; i++)
   {  for (j=0; j<G.row_number[i]; j++)
      {  u[i] += G.value[i][j] * v[G.index[i][j]];
         }
      }
   dwtLin(u,F,M,p,np);
   for (i=0; i<np; i++) rhs[i] += 4*pi*u[i]/(epsilon-1);
   memset(u,0,np*sizeof(double));
   i = WEMPCG_pwl(&S_i,rhs,u,eps,W,F,p,M);
   printf("Solving the 2nd linear system:   %d iterations\n",i);
   }
else if ((CASE == 3) || (CASE == 4))
{  WEMRHS_pwl2(&rhs,W,E,T,p,M,np);
   i = WEMPCG_pwl(&S_i,rhs,u,eps,W,F,p,M);		/* u = V_i^(-1)*N_f */
   printf("Solving the 1st linear system:   %d iterations\n",i);
   memset(rhs,0,np*sizeof(double));
   for (i=0; i<np; i++)			   	/* rhs = V_e*u */
   {  for (j=0; j<S_e.row_number[i]; j++) rhs[i] += S_e.value1[i][j] * u[S_e.index[i][j]];
      }
   i = WEMPGMRES_pwl3(&S_i,&S_e,rhs,v,eps,W,F,p,M);	/* solve complicated_system u = A^(-1)*rhs */ 
   printf("Solving the 2nd linear system:   %d iterations\n",i);   
   for (i=0; i<np; i++) u[i] -= 4*pi*v[i]; 	/* u = u - 4*pi*v */ 
   }

time(&t2);
printf("Computation time:                %g secs.\n\n",difftime(t2,t3));

/* Energie-Berechnung */
tdwtLin(u,F,M,p,np);
res = energy(u,F,T,p,M);
time(&t2);
printf("Over-all computation time:       %g secs.\n",difftime(t2,t1));

/* Speicherplatz freigeben */
printf("\n\n\n");
if ((CASE == 3) || (CASE == 4)) free_sparse2(&S_e);
free_waveletlist(&W,np);
free_elementlist(&E,p,M);
free_interpolate(&T,p,M);
free_patchlist(&F,nf);
free_sparse2(&S_i);
free(rhs);
free(P);
free(u);
free(v);
return(0);
}
