/**************
 *  WEMPCG_pwl.c  *
 **************/


/*=============================================================*
 *  WEMPCG_pwl(A,b,x,epsi,W,F,p,M)                                 *
 *	                                                       *
 *  Verfahren der konjugierten Gradienten zur Loesung des      *
 *  linearen Gleichungssystems                                 *
 *	                                                       *
 *	      		A1*x = b.               	       *
 *	                                                       *
 *  Vorkonditionierung per Wavelet-Preconditioner.             *
 *	                                                       *
 *  Parameter :                                                *
 *		A    : Matrix im sparse2-Format                *
 *		G    : Massenmatrix                            *
 *		b    : rechte Seite                            *
 *		x    : Startwert und Endwert                   *
 *		epsi : Genauigkeit                             *
 *		W    : Liste der Wavelets                      *
 *		F    : Elementliste der Einskalenbasis         *
 *		p    : Zahl der Paramtergebiete                *
 *		M    : Zahl der Level                          *
 *=============================================================*/


#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "sparse.h"
#include "sparse2.h"
#include "intvector_pwl.h"
#include "vector3.h"
#include "basis.h"
#include "kern.h"
#include "dwt_pwl.h"
#include "precond.h"
#include "WEMPCG_pwl.h"


unsigned int WEMPCG_pwl(A,b,x,epsi,W,F,p,M)
sparse2 	*A;
double          *b, *x, epsi;
wavelet		*W;
unsigned int	**F, p, M;
{
unsigned int    i, j, k;
double  	u, v, omg;
double  	*z, *r, *d, *Ad;
sparse		G;

/* Speicherplatz allokieren */
r  = (double*) calloc(A->n,sizeof(double));
d  = (double*) calloc(A->n,sizeof(double));
z  = (double*) malloc(A->n*sizeof(double));
Ad = (double*) malloc(A->n*sizeof(double));

/* berechne Gram'sche Matrix */
init_sparse(&G,A->n,A->n,10);
single_scale_gram(&G,F,p,M);

/* r = b-A1*x */
memcpy(r,b,A->n*sizeof(double));
for (i=0; i<A->n; i++)
{  for (j=0; j<A->row_number[i]; j++)
   {  r[i] -= A->value1[i][j] * x[A->index[i][j]];
      }
   }

/* d = precond(r,M) */
precond(d,r,&G,W,F,p,M);

/* u = (r,d) */
u = 0;
for (i=0; i<A->n; i++) u += r[i] * d[i];

/* Iteration */
for (k=0; sqrt(u)>epsi; k++)
{  
   /* Ad = A*d */
   memset(Ad,0,A->n*sizeof(double));
   for (i=0; i<A->n; i++)
   {  for (j=0; j<A->row_number[i]; j++)
      {  Ad[i] += A->value1[i][j] * d[A->index[i][j]];
         }
      }

   /* omg = u / (d,Ad) und v = u */
   omg = 0;
   for (i=0; i<A->n; i++) omg += d[i] * Ad[i];
   omg = u/omg;
   v = u;

   /* x = x + omg * d, r = r - omg * Ad */
   for (i=0; i<A->n; i++)
   {  x[i] += omg *  d[i];
      r[i] -= omg * Ad[i];
      }

   /* z = precond(r,M) */
   precond(z,r,&G,W,F,p,M);

   /* u = (r,z) */
   u = 0;
   for (i=0; i<A->n; i++) u += r[i] * z[i];

   /* d = z + u/v * d */
   omg = u/v;
   for (i=0; i<A->n; i++) d[i] = z[i] + omg * d[i];
   }

free_sparse(&G);
free(Ad);
free(r);
free(d);
free(z);
return(k);
}
