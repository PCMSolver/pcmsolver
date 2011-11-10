/****************
 *  PostProc.c  *
 ****************/


/*===============================================*
 *  Postprocessing: setzt alle Eintraege, deren  *
 *  Betrag kleiner als ein levelabhaengiger      *
 *  Abschneideparameter ist, auf Null.           *
 *===============================================*/


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "constants.h"
#include "sparse2.h"
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "postproc.h"


double postproc(T,W,E,p,M)
sparse2         *T;	      	/* komprimierte Steifigkeitsmatrix */
wavelet		*W;		/* Liste der Wavelets              */
element		*E;		/* hierarchische Elementliste      */
unsigned int	p;	      	/* Anzahl der Patches              */
unsigned int	M;	      	/* 2^M*2^M Elemente pro Patch      */
{
signed int    	i, j, k;	/* Laufindizes                     */
double		**c;		/* Abschneideparameter 		   */
double		*D1, *D2;	/* Diagonaleintrage von T          */
double		*v1, *v2;	/* neue Wertevektoren              */
unsigned int	*w;		/* neuer Indexvktor                */
unsigned int    nnz;	      	/* nunber of nonzero elements of T */

/* besimme die Diagonalen von T */
D1 = (double*) malloc(T->n*sizeof(double));
D2 = (double*) malloc(T->n*sizeof(double));
for (i=0; i<T->n; i++) 
{  get_sparse2(T,i,i,&D1[i],&D2[i]);
   D1[i] = sqrt(fabs(D1[i]));
   D2[i] = sqrt(fabs(D2[i]));
   }

/* berechne Abschneideparameter */
c = (double**) malloc((M+1)*sizeof(double*));
for (i=0; i<=M; i++)
{  c[i] = (double*) malloc((M+1)*sizeof(double));
   for (j=0; j<=M; j++)
   {  c[i][j] = pow(0.5,(2*M-(i+j))*(2*dp-op)/(2*td+op));
      if (c[i][j] > pow(0.5,fabs(i-j))) c[i][j] = pow(0.5,fabs(i-j));
      c[i][j] *= b * pow(0.5,dp*(2*M-(i+j)));
      }
   }

/* a-posteriori-Kompression */
nnz = 0;
for (i=0; i<T->n; i++)
{  
   /* Suche relavante Matrixeintraege aus */
   k=0;
   v1 = (double*) malloc(T->row_number[i]*sizeof(double));
   v2 = (double*) malloc(T->row_number[i]*sizeof(double));
   w = (unsigned int*) malloc((T->row_number[i]+1)*sizeof(unsigned int));

   for (j=0; j<T->row_number[i]; j++)
   {  if ( (fabs(T->value1[i][j]) >= D1[i]*D1[T->index[i][j]]*c[W[i].level][W[T->index[i][j]].level]) ||
           (fabs(T->value2[i][j]) >= D2[i]*D2[T->index[i][j]]*c[W[i].level][W[T->index[i][j]].level]) )
      {  v1[k] = T->value1[i][j];
         v2[k] = T->value2[i][j];
	 w[k]  = T->index[i][j];
	 k++;
	 }
      }

   /* Eintraege kopieren und Speicherplatz freigeben */
   free(T->index[i]);
   free(T->value1[i]);
   free(T->value2[i]);
   T->max_row_number[i] = T->row_number[i] = k;
   T->value1[i] = (double*)       realloc(v1, k   *sizeof(double));
   T->value2[i] = (double*)       realloc(v2, k   *sizeof(double));
   T->index[i] = (unsigned int*) realloc(w,(k+1)*sizeof(unsigned int));
   T->index[i][k] = T->n;
   nnz += k;
   }
   
/* Speicherplatz wieder freigeben */
/*printf("A-posteriori compression:        %.5f %%\n",100.0*nnz/T->n/T->n);*/
for (i=0; i<=M; i++) free(c[i]);
free(D1);
free(D2);
free(c);
return 100.0*nnz/T->n/T->n;
}
