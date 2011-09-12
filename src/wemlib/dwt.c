/***********
 *  dwt.c  *
 ***********/
 
 
/*===========================================*
 *  Dieses Modul enthaelt alle Routinen der  *
 *  schnellen Wavelet-Transformationen.      *
 *===========================================*/
 
 
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "sparse.h"
#include "mask.h"
#include "dwt.h"


void dwtKon(a,M,nf)
/* Diskrete Wavelet-Transformation */
double		*a;	        /* gesuchter/gegebener Koeffizientenvektor */
unsigned int	M;	        /* (2^M*2^M) Elemente pro Patch            */
unsigned int	nf;	        /* number of elements  		           */
{
unsigned int 	p;		/* Anzhl der Parametergebiete		   */
unsigned int	m;		/* Laufindex fuer Level			   */
unsigned int	n;		/* Elemente pro Patch auf Level m          */
unsigned int	i1, i2, i3;     /* Laufindizes fuer die Elemente           */
unsigned int	zg;		/* Nullelement von Patch i1 auf Level m    */
sparse		T, L;		/* Maskenmatrizen                          */
unsigned int	s, t;		/* Laufindizes durch die Maskenmatrizen    */
double		*b;		/* Hilfsvektor 			           */

/* SCHLEIFE UEBER DIE GITTER */
n = 1 << M;
p = nf/(n*n);
b = (double*) calloc(nf,sizeof(double));
for (m=M; m>=minLevel; m--)
{  
   dwt_mask(&T,&L,m);	/* berechne die Masken T und L  */
   n = 1 << (m-1);	/* p*n*n Elemente auf Level m-1 */
   
   /* 1. bilde Skalierungsfunktionen und Wavelets */
   for (i1=0; i1<p; i1++)
   {  zg = i1*n*n;	/* Nullpunkt von Patch i1 auf Level m-1 */
      for (i2=0; i2<n; i2++)
      {  for (i3=0; i3<n; i3++)
         {  
	    /* Skalierungsfunktion */
	    b[zg+n*i2+i3] = 0.5 * ( a[4*zg+(4*n*i2    )+(2*i3  )] \
	                          + a[4*zg+(4*n*i2    )+(2*i3+1)] \
			          + a[4*zg+(4*n*i2+2*n)+(2*i3  )] \
			          + a[4*zg+(4*n*i2+2*n)+(2*i3+1)]);

	    /* Wavelet: psi(s)*phi(t) */
	    for (s=0; s<T.row_number[i3]; s++)
	    {  b[zg+p*n*n+n*i2+i3] += 0.5 * T.value[i3][s] * ( a[4*zg+(4*n*i2    )+T.index[i3][s]] \
	                                                     + a[4*zg+(4*n*i2+2*n)+T.index[i3][s]] );
               }

	    /* Wavelet: phi_1(s)*psi(t) und phi_2(s)*psi(t) */
	    for (t=0; t<T.row_number[i2]; t++)
	    {  b[zg+2*p*n*n+n*i2+i3] += 0.5 * T.value[i2][t] * a[4*zg+(2*n*T.index[i2][t])+(2*i3  )];
               b[zg+3*p*n*n+n*i2+i3] += 0.5 * T.value[i2][t] * a[4*zg+(2*n*T.index[i2][t])+(2*i3+1)];
               }
       	    }
	 }
      }

   /* 2. Skalierungsfunktionen nach a kopieren */
   memcpy(a,b,p*n*n*sizeof(double));
   memset(b,0,p*n*n*sizeof(double));
   free_sparse(&T);
   }

/* Vektor b wieder freigeben */
memcpy(&a[p*n*n],&b[p*n*n],(nf-p*n*n)*sizeof(double));
free(b);
return;
}


void tdwtKon(a,M,nf)
/* transformierte Diskrete Wavelet-Transformation */
double		*a;	        /* gesuchter/gegebener Koeffizientenvektor */
unsigned int	M;	        /* (2^M*2^M) Elemente pro Patch            */
unsigned int	nf;	        /* number of elements  		           */
{
unsigned int 	p;		/* Anzhl der Parametergebiete		   */
unsigned int	m;		/* Laufindex fuer Level			   */
unsigned int	n;		/* Elemente pro Patch auf Level m          */
unsigned int	i1, i2, i3;     /* Laufindizes fuer die Elemente           */
unsigned int	zg;		/* Nullelement von Patch i1 auf Level m    */
sparse		T, L;		/* Maskenmatrizen                          */
unsigned int	s, t;		/* Laufindizes durch die Maskenmatrizen    */
double		*b;		/* Hilfsvektor 			           */

/* SCHLEIFE UEBER DIE GITTER */
p = nf/(1<<2*M);
b = (double*) calloc(nf,sizeof(double));
for (m=minLevel; m<=M; m++)
{
   dwt_mask(&T,&L,m);	/* berechne die Masken T und L  */
   n = 1 << (m-1);	/* p*n*n Elemente auf Level m-1 */
   
   /* 1. Skalierungsfunktionen und Wavelets bilden */
   for (i1=0; i1<p; i1++)
   {  zg = i1*n*n;	/* Nullpunkt von Patch i1 auf Level m */
      for (i2=0; i2<n; i2++)
      {  for (i3=0; i3<n; i3++)
         {  
	    b[4*zg+(4*n*i2    )+2*i3  ] += 0.5 * a[zg+n*i2+i3];
	    b[4*zg+(4*n*i2    )+2*i3+1] += 0.5 * a[zg+n*i2+i3];
	    b[4*zg+(4*n*i2+2*n)+2*i3  ] += 0.5 * a[zg+n*i2+i3];
	    b[4*zg+(4*n*i2+2*n)+2*i3+1] += 0.5 * a[zg+n*i2+i3];
	    
	    for (s=0; s<T.row_number[i3]; s++)
	    {  b[4*zg+(4*n*i2    )+T.index[i3][s]] += 0.5 * T.value[i3][s] * a[zg+p*n*n+n*i2+i3];
	       b[4*zg+(4*n*i2+2*n)+T.index[i3][s]] += 0.5 * T.value[i3][s] * a[zg+p*n*n+n*i2+i3];
	       }
	       
	    for (t=0; t<T.row_number[i2]; t++)
	    {  b[4*zg+(2*n*T.index[i2][t])+(2*i3  )] += 0.5 * T.value[i2][t] * a[zg+2*p*n*n+n*i2+i3];
	       b[4*zg+(2*n*T.index[i2][t])+(2*i3+1)] += 0.5 * T.value[i2][t] * a[zg+3*p*n*n+n*i2+i3];
	       }
       	    }
	 }
      }

   /* 2. Skalierungsfunktionen nach a kopieren */
   memcpy(a,b,4*p*n*n*sizeof(double));
   memset(b,0,4*p*n*n*sizeof(double));
   free_sparse(&T);
   }

/* Vektor b wieder freigeben */
free(b);
return;
}
