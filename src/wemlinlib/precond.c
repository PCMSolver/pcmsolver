/*******************
 *  precond_mod.c  *
 *******************/
 

/*========================================================*
 *  Berechnet den Preconditioner bzgl. des modifizierten  *
 *  Skalarproduktes und definiert dessen Anwendung auf    *
 *  einen Vektor.                                         *
 *========================================================*/
 
 
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "sparse.h"
#include "basis.h"  
#include "precond.h"
#include "dwt.h"


void inv_A_times_x(sparse *A, double *x, unsigned int **F, unsigned int p, unsigned int M);


void single_scale_gram(G,F,p,M)
/* berechnet die Massenmatrix in der Einskalenbasis. */
sparse		*G;             /* zu berechnende Gram-Matrix im sparse-Format */
unsigned int	**F;		/* Elementliste der Einskalenbasis             */
unsigned int    p;		/* Anzahl der Patches		               */
unsigned int	M;		/* (2^M*2^M) Elemente pro Patch                */
{
unsigned int	N = 1 << M;     /* N*N Elemente pro Patch                      */
unsigned int	i, j;		/* Laufindizes		                       */
double		c[16];		/* Integralgroessen 			       */

/* Initialisierung */
c[0] = c[5] = c[10] = c[15] = 1./9;
c[1] = c[4] = c[ 9] = c[12] = 1./18;
c[2] = c[7] = c[ 8] = c[13] = 1./36;
c[3] = c[6] = c[11] = c[14] = 1./18;

/* Aufbau der Massenmatrix */
for (i=0; i<p*N*N; i++)
{  for (j=0; j<16; j++) add_sparse(G,F[i][j%4],F[i][j/4],c[j]);	/* L^2-normiert! */
   }

return;
}


void inv_A_times_x(A,x,F,p,M)
/* Loest das lineare Gleichungssystem T*A*T'*y = x,
   wobei der Nullvektor als Startvektor verwendet 
   und die Loesung y in x geschrieben wird. */
sparse 		*A;
double          *x;
unsigned int	**F;
unsigned int	p, M;
{
unsigned int    i, j;
double  	u, u_start, v, omg;
double  	*r, *d, *z, *Ad;

/* Speicherplatz allokieren */
z  = (double*) malloc(A->n*sizeof(double));
r  = (double*) calloc(A->n,sizeof(double));
d  = (double*) malloc(A->n*sizeof(double));
Ad = (double*) malloc(A->n*sizeof(double));

/* r = x - T*A*T'*x, d = r und u = (r,r) */
u = 0;
memcpy(d,x,A->n*sizeof(double));
tdwtLin(d,F,M,p,A->n);
for (i=0; i<A->n; i++)
{  for (j=0; j<A->row_number[i]; j++) r[i] -= A->value[i][j] * d[A->index[i][j]];
   }
dwtLin(r,F,M,p,A->n);
for (i=0; i<A->n; i++)
{  r[i] += x[i];
   d[i] = r[i];
   u += r[i] * r[i];
   }
u_start = u;

/* Iteration */
while (sqrt(u/u_start)>eps)
{
   /* Ad = A*d */
   memcpy(z,d,A->n*sizeof(double));
   memset(Ad,0,A->n*sizeof(double));
   tdwtLin(z,F,M,p,A->n);
   for (i=0; i<A->n; i++)
   {  for (j=0; j<A->row_number[i]; j++) Ad[i] += A->value[i][j] * z[A->index[i][j]];
      }
   dwtLin(Ad,F,M,p,A->n);

   /* omg = u / (d,Ad) und v = u */
   omg = 0;
   for (i=0; i<A->n; i++) omg += d[i] * Ad[i];
   omg = u/omg;
   v = u;

   /* x = x + omg * d, r = r - omg * Ad und u = (r,r) */
   u = 0;
   for (i=0; i<A->n; i++)
   {  x[i] += omg * d[i];
      r[i] -= omg * Ad[i];
      u += r[i] * r[i];
      }

   /* d = r +  u/v * d */
   omg = u/v;
   for (i=0; i<A->n; i++) d[i] = r[i] + omg * d[i];
   }

free(Ad);
free(z);
free(r);
free(d);
return;
}


void precond(a,b,G,W,F,p,M)
/* Anwendung des Preconditioners auf den Vektor b, der 
   NICHT veraendert wird. Das Ergebnis wird in a gespeichert */
double		*a;		/* gesuchter Koeffizientenvektor   */
double		*b;		/* gegebener Koeffizientenvektor   */
sparse		*G;		/* Gram'sche Matrix                */
wavelet		*W;		/* Liste der Wavelets              */
unsigned int	**F;		/* Elementliste der Einskalenbasis */
unsigned int	p;		/* Anzahl der Patches              */
unsigned int	M;		/* (2^M*2^M) Elemente pro Patch    */
{
unsigned int	np = G->n;	/* number of basis functions	   */
unsigned int	i;		/* Laufindex durch ein Level       */

/* wende Wavelet-Preconditioner an */
for (i=0; i<np; i++) a[i] = pow(2,0.5*W[i].level) * b[i];

/* Multipliziere mit der Wavelet-Massenmatrix */
inv_A_times_x(G,a,F,p,M);

/* wende Wavelet-Preconditioner an */
for (i=0; i<np; i++) a[i] *= pow(2,0.5*W[i].level);

return;
}
