/***************
 *  precond.c  *
 ***************/
 

/*========================================================*
 *  Berechnet den Preconditioner bzgl. des modifizierten  *
 *  Skalarproduktes und definiert dessen Anwendung auf    *
 *  einen Vektor.                                         *
 *========================================================*/
 
 
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "vector2.h"
#include "vector3.h" 
#include "dwt.h"
#include "precond.h"


void inv_G_times_x(x,p,M)
/* Loest das lineare Gleichungssystem T*A*T'*y = x, wobei x als 
   Startvektor verwendet und die Loesung y in x geschrieben wird. */
double          *x;
unsigned int    p;
unsigned int	M;
{
unsigned int    i, nf;
double  	u, u_start, v, omg;
double  	*r, *d, *Ad;

/* Speicherplatz allokieren */
nf = p*(1<<M)*(1<<M);
r  = (double*) calloc(nf,sizeof(double));
d  = (double*) malloc(nf*sizeof(double));
Ad = (double*) malloc(nf*sizeof(double));

/* r = x - T*T'*x, d = r und u = (r,r) */
u = 0;
for (i=0; i<nf; i++) r[i] = -x[i];
tdwtKon(r,M,nf);
dwtKon(r,M,nf);
for (i=0; i<nf; i++)
{  r[i] += x[i];
   d[i] = r[i];
   u += r[i] * r[i];
   }

 if(sqrt(u) > eps) {
     u_start = u;

/* Iteration */
while (sqrt(u/u_start)>eps)
{
   /* Ad = A*d */
   memcpy(Ad,d,nf*sizeof(double));
   tdwtKon(Ad,M,nf);
   dwtKon(Ad,M,nf);

   /* omg = u / (d,Ad) und v = u */
   omg = 0;
   for (i=0; i<nf; i++) omg += d[i] * Ad[i];
   omg = u/omg;
   v = u;

   /* x = x + omg * d, r = r - omg * Ad und u = (r,r) */
   u = 0;
   for (i=0; i<nf; i++)
   {  x[i] += omg * d[i];
      r[i] -= omg * Ad[i];
      u += r[i] * r[i];
      }

   /* d = r +  u/v * d */
   omg = u/v;
   for (i=0; i<nf; i++) d[i] = r[i] + omg * d[i];
 }
 } else {
     for(i=0; i < nf; i++) {
         x[i] = 0.0l;
     }
 }

free(Ad);
free(r);
free(d);
return;
}


void precond(a,b,p,M)
/* Anwendung des Preconditioners auf den Vektor b, der 
   NICHT veraendert wird. Das Ergebnis wird in a gespeichert */
double		*a;		/* gesuchter Koeffizientenvektor  */
double		*b;		/* gegebener Koeffizientenvektor  */
unsigned int    p;		/* Anzahl der Patches		  */
unsigned int	M;		/* (2^M*2^M) Elemente pro Patch   */
{
unsigned int	m;		/* Laufindex durch die Level      */
unsigned int	min_i;		/* kleinster Index im Level m     */
unsigned int	max_i;		/* groesster Index im Level m     */
unsigned int	i;		/* Laufindex durch ein Level      */
double		c;		/* Wichtungsfaktor eines Wavelets */

/* wende Wavelet-Preconditioner an */
max_i = 0;
for (m=minLevel; m<=M; m++)
{  min_i = max_i;
   max_i = p*(1<<m)*(1<<m);
   c = pow(2,-0.5*op*m);
   for (i=min_i; i<max_i; i++) a[i] = c * b[i];
   }

/* Multipliziere mit der Wavelet-Massenmatrix */
inv_G_times_x(a,p,M);

/* wende Wavelet-Preconditioner an */
max_i = 0;
for (m=minLevel; m<=M; m++)
{  min_i = max_i;
   max_i = p*(1<<m)*(1<<m);
   c = pow(2,-0.5*op*m);
   for (i=min_i; i<max_i; i++) a[i] *= c;
   }
return;
}
