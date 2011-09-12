/*******************
 *  Compression.c  *
 *******************/


/*===========================================*
 *  Bestimmt in der Steifigkeitsmatrix alle  *
 *  nach der 1. und 2. Kompression zu        *
 *  berechnenden Eintraege.                  *
 *===========================================*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "sparse2.h"
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "compression.h"


/*========================================================*
 *  Deklaration der Hilfsfunktionen und der bounding box  *
 *========================================================*/

typedef struct		/* Typdeklaration bounding box */
{  double	r_x, r_y, r_z;
   double	m_x, m_y, m_z;
   }
bounding_box_gamma;


typedef struct		/* Typdeklaration bounding box */
{  unsigned int	patch;
   double	r_x, r_y;
   double	m_x, m_y;
   }
bounding_box_square;


void compute_bounding_boxes_gamma(bounding_box_gamma **B, wavelet *W, element *E, unsigned int p, unsigned int M);


void compute_bounding_boxes_square(bounding_box_square **B, wavelet *W, element *E, unsigned int p, unsigned int M);


unsigned int wavelet_wavelet_criterion(bounding_box_gamma *B1, bounding_box_square *B2,
	wavelet *W, element *E, unsigned int ind1, unsigned int ind2, double c1, double c2);


/*==================================================*
 *  Berechnet die bounding boxes fuer die Wavelets  *
 *==================================================*/

void compute_bounding_boxes_gamma(B,W,E,p,M)
bounding_box_gamma	**B;			/* bounding boxes             */
wavelet			*W;			/* Liste der Wavelets         */
element         	*E;			/* hierarchische Elementliste */
unsigned int		p;			/* Anzahl der Patches         */
unsigned int		M;			/* Zahl der Level             */
{
unsigned int		nf;			/* Anzahl der Basisfunktionen */
unsigned int		i, j, k;		/* Laufindizes                */
double			min_x, min_y, min_z;	/* minimale Ausdehnung        */
double			max_x, max_y, max_z;	/* maximale Ausdehnung        */

nf = p*(1<<2*M);
(*B) = (bounding_box_gamma*) malloc(nf*sizeof(bounding_box_gamma));
for (i=0; i<nf; i++)
{  k = W[i].element[0];
   min_x = E[k].midpoint.x-E[k].radius;
   max_x = E[k].midpoint.x+E[k].radius;
   min_y = E[k].midpoint.y-E[k].radius;
   max_y = E[k].midpoint.y+E[k].radius;
   min_z = E[k].midpoint.z-E[k].radius;
   max_z = E[k].midpoint.z+E[k].radius;
   for (j=1; j<W[i].element_number; j++)
   {  k = W[i].element[j];
      if (E[k].midpoint.x-E[k].radius < min_x) min_x = E[k].midpoint.x-E[k].radius;
      if (E[k].midpoint.y-E[k].radius < min_y) min_y = E[k].midpoint.y-E[k].radius;
      if (E[k].midpoint.z-E[k].radius < min_z) min_z = E[k].midpoint.z-E[k].radius;
      if (E[k].midpoint.x+E[k].radius > max_x) max_x = E[k].midpoint.x+E[k].radius;
      if (E[k].midpoint.y+E[k].radius > max_y) max_y = E[k].midpoint.y+E[k].radius;
      if (E[k].midpoint.z+E[k].radius > max_z) max_z = E[k].midpoint.z+E[k].radius;
      }
   (*B)[i].r_x = 0.5*(max_x-min_x);
   (*B)[i].r_y = 0.5*(max_y-min_y);
   (*B)[i].r_z = 0.5*(max_z-min_z);
   (*B)[i].m_x = 0.5*(max_x+min_x);
   (*B)[i].m_y = 0.5*(max_y+min_y);
   (*B)[i].m_z = 0.5*(max_z+min_z);
   }
return;
}   


void compute_bounding_boxes_square(B,W,E,p,M)
bounding_box_square	**B;		/* bounding boxes             */
wavelet			*W;		/* Liste der Wavelets         */
element         	*E;		/* hierarchische Elementliste */
unsigned int		p;		/* Anzahl der Patches         */
unsigned int		M;		/* Zahl der Level             */
{
unsigned int		nf;		/* Anzahl der Basisfunktionen */
unsigned int		i, j, k, l;	/* Laufindizes                */
double			min_x, min_y;	/* minimale Ausdehnung        */
double			max_x, max_y;	/* maximale Ausdehnung        */

l = p;
nf = p*(1<<2*M);
(*B) = (bounding_box_square*) malloc(nf*sizeof(bounding_box_square));
for (i=0; i<nf; i++)
{  k = W[i].element[0];
   (*B)[i].patch = E[k].patch;
   min_x =  E[k].index_s   *1./(1<<E[k].level);
   max_x = (E[k].index_s+1)*1./(1<<E[k].level);
   min_y =  E[k].index_t   *1./(1<<E[k].level);
   max_y = (E[k].index_t+1)*1./(1<<E[k].level);
   for (j=1; j<W[i].element_number; j++)
   {  k = W[i].element[j];
      if (E[k].patch != (*B)[i].patch) 
      {  (*B)[i].patch = l++;
         break;		/* Wavelet liegt auf mehreren Patches */
	 }
      if (E[k].index_s   < min_x*(1<<E[k].level)) min_x =  E[k].index_s   *1./(1<<E[k].level);
      if (E[k].index_t   < min_y*(1<<E[k].level)) min_y =  E[k].index_t   *1./(1<<E[k].level);
      if (E[k].index_s+1 > max_x*(1<<E[k].level)) max_x = (E[k].index_s+1)*1./(1<<E[k].level);
      if (E[k].index_t+1 > max_y*(1<<E[k].level)) max_y = (E[k].index_t+1)*1./(1<<E[k].level);
      }
   (*B)[i].r_x = 0.5*(max_x-min_x);
   (*B)[i].r_y = 0.5*(max_y-min_y);
   (*B)[i].m_x = 0.5*(max_x+min_x);
   (*B)[i].m_y = 0.5*(max_y+min_y);
   }
return;
}   


/*=========================================================*
 *  wavelet_wavelet_criterion: 1.+2. Kompression	   *
 *  result = 0: T(ind,ind2) = 0 gemaess 1./2. Kompression  *
 *  result = 1: sonst 					   *
 *=========================================================*/

unsigned int wavelet_wavelet_criterion(B1,B2,W,E,ind1,ind2,c1,c2)
bounding_box_gamma	*B1;		/* bounding boxes                                */
bounding_box_square	*B2;		/* bounding boxes                                */
wavelet			*W;		/* Liste der Wavelets                            */
element         	*E;		/* hierarchische Elementliste                    */
unsigned int		ind1, ind2;	/* Indizes der beiden Wavelets                   */
double			c1, c2;		/* Abschneideparameter 1./2.  Kompression        */
{
unsigned int		i, j;		/* Laufindizes durch die Elemente                */
unsigned int		k, l;		/* Indizes der Elemente                          */
double			h1, h2;		/* halbe Kantenlaenge von Element1 bzw. Element2 */
double			s1, t1;		/* (s1,t1) = Mittelpunkt von Element1            */
double			s2, t2;		/* (s2,t2) = Mittelpunkt von Element2            */
double			dist;		/* Abstand zwischen den MP zweier Elemente       */
double			dx, dy, dz;	/* x/y/z-Abstand zweier Wavelets                 */

/* 1. Kompression: nehme die bounding boxes um die Wavelets */
if (ind1 < ind2) return(0);
dx = fabs(B1[ind1].m_x-B1[ind2].m_x)-B1[ind1].r_x-B1[ind2].r_x;
dy = fabs(B1[ind1].m_y-B1[ind2].m_y)-B1[ind1].r_y-B1[ind2].r_y;
dz = fabs(B1[ind1].m_z-B1[ind2].m_z)-B1[ind1].r_z-B1[ind2].r_z;
if (dx < 0) dx = 0;
if (dy < 0) dy = 0;
if (dz < 0) dz = 0;
if (dx*dx+dy*dy+dz*dz >= c1*c1) return(0);

/* 2. Kompression: beide Wavelets leben auf dem selben Patch */
if (B2[ind1].patch == B2[ind2].patch)
{  dx = fabs(B2[ind1].m_x-B2[ind2].m_x)-B2[ind1].r_x-B2[ind2].r_x;
   dy = fabs(B2[ind1].m_y-B2[ind2].m_y)-B2[ind1].r_y-B2[ind2].r_y;
   if (dx < 0) dx = 0;
   if (dy < 0) dy = 0;
   if (dx+dy > 0) 	/* Wavelets sind disjunkt */
   {  if (dx*dx+dy*dy >= c2*c2) return(0); else return(1);
      for (i=0; i<W[ind1].element_number; i++)
      {  for (j=0; j<W[ind2].element_number; j++)
         {  if (distance(&E[W[ind1].element[i]],&E[W[ind2].element[j]]) < c2) return(1);
            }
         }
      }
   else		/* Wavelets sind nicht disjunkt */
   {  /* teste, ob W[ind1]+c2 \not\subset W[ind2] */
      if  ( (B2[ind1].m_x-B2[ind1].r_x-c2 < B2[ind2].m_x-B2[ind2].r_x) \
         || (B2[ind1].m_x+B2[ind1].r_x+c2 > B2[ind2].m_x+B2[ind2].r_x) \
         || (B2[ind1].m_y-B2[ind1].r_y-c2 < B2[ind2].m_y-B2[ind2].r_y) \
         || (B2[ind1].m_y+B2[ind1].r_y+c2 > B2[ind2].m_y+B2[ind2].r_y) ) return(1);

      /* teste bounding box von W[ind1] gegen Elemente von W[ind2] */
      for (j=0; j<W[ind2].element_number; j++)
      {  l = W[ind2].element[j];
         h2 = 0.5/(1<<E[l].level);
         s2 = h2*(2*E[l].index_s + 1);
         t2 = h2*(2*E[l].index_t + 1);
         
	 dx = fabs(B2[ind1].m_x-s2)-B2[ind1].r_x-h2;
         dy = fabs(B2[ind1].m_y-t2)-B2[ind1].r_y-h2;
         if (dx < 0) dx = 0;
         if (dy < 0) dy = 0;
         if (dx+dy > 0) 	/* Wavelet und Element sind disjunkt */
         {  if (dx*dx+dy*dy < c2*c2) return(1);
	    }
	 else			/* Wavelet und Element sind nicht disjunkt */   
         {  if  ( (B2[ind1].m_x-B2[ind1].r_x-c2 < s2-h2) \
               || (B2[ind1].m_x+B2[ind1].r_x+c2 > s2+h2) \
               || (B2[ind1].m_y-B2[ind1].r_y-c2 < t2-h2) \
               || (B2[ind1].m_y+B2[ind1].r_y+c2 > t2+h2) ) return(1);
            }
         }
      }
   return(0);
   }

/* 2. Kompression: Wavelets nicht auf demseleben Patch -> teste einzelne Elemente */
for (i=0; i<W[ind1].element_number; i++)
{  k = W[ind1].element[i];
   h1 = 0.5/(1<<E[k].level);
   s1 = h1*(2*E[k].index_s + 1);
   t1 = h1*(2*E[k].index_t + 1);

   for (j=0; j<W[ind2].element_number; j++)
   {  l = W[ind2].element[j];
   
      /* entweder Elemente liegen auf dem selben Patch */
      if (E[k].patch == E[l].patch)
      {  h2 = 0.5/(1<<E[l].level);
         s2 = h2*(2*E[l].index_s + 1);
         t2 = h2*(2*E[l].index_t + 1);

	 /* berechne Unendlichnorm des Abstandes der beiden Mittelpunkte */
	 dist = (fabs(s1-s2) < fabs(t1-t2)) ? fabs(t1-t2) : fabs(s1-s2);
         
	 /* Abstand der Elemente */
	 if (fabs(fabs(dist-h2)-h1) < c2) return(1);
         }

      /* oder Elemente liegen auf verschiedenen Patches */
      else if (distance(&E[k],&E[l]) < c2) return(1);
      }
   }
return(0);
}


/*=================*
 *  Hauptprogramm  *
 *=================*/
  
double compression(T,W,E,p,M)
sparse2          	*T;	      	/* komprimierte Steifigkeitsmatrix               */
wavelet			*W;		/* Liste der Wavelets                            */
element			*E;		/* hierarchische Elementliste                    */
unsigned int		p;		/* Zahl der Patches                              */
unsigned int		M;	      	/* 2^M*2^M Patches pro Parametergebiet           */
{
unsigned int		nf;		/* Anzahl der Basisfunktionen                    */
bounding_box_gamma	*B1;		/* bounding boxes fuer die Wavelets              */
bounding_box_square	*B2;		/* bounding boxes fuer die Wavelets              */
unsigned int 		m1, m2;		/* Laufindizes durch Level		         */
unsigned int		ind1, ind2;	/* Argumente der Vater-Wavelets                  */
unsigned int		i, j;		/* Spaltenlaufindizes durch Steifigkeitsmatrix T */
unsigned int		k, l;		/* Laufindizes durch Soehne                      */
unsigned int		*rn;		/* Zaehlvektor zur Symmetrisierung		 */
double			max_radius;	/* maximaler Umkreis der Elemente am Level 0	 */
double			**c1, **c2;	/* Abschneideparameter in 1./2. Kompression      */
double			d1, d2;		/* Hilfskonstanten                               */
unsigned int    	nnz;	      	/* nunber of nonzero elements of T      	 */

/* berechne Abschneideparameter */
max_radius = 0;
for (i=0; E[i].level==0; i++) if (max_radius < E[i].radius) max_radius = E[i].radius;
c1 = (double**) malloc((M+1)*sizeof(double*));
c2 = (double**) malloc((M+1)*sizeof(double*));
for (i=0; i<=M; i++)
{  c1[i] = (double*) malloc((i+1)*sizeof(double));
   c2[i] = (double*) malloc((i+1)*sizeof(double));
   for (j=0; j<=i; j++)
   {  c1[i][j] = a*pow(2,(M*(2*dp-op)-(i+j)*(dp+td))/(2*td+op));
      d1 = pow(2,(M*(2*dp-op)-(i+j)*dp-i*td)/(td+op));		/* alter Parameter */
      d2 = pow(2,(2*M-(i+j))*(2*dp-op)/((2*td+op)*(td+op)))*pow(2,(M*(2*dp-op)-i*(dp+td+1)-j*(dp-1))/(td+op));
      if (d1 > d2) c2[i][j] = a*d1; else c2[i][j] = a*d2;
      if (c1[i][j] < a/(1<<j)) c1[i][j] = a/(1<<j);
      if (c2[i][j] < a/(1<<i)) c2[i][j] = a/(1<<i);
      if (c1[i][j] > c2[i][j]) c1[i][j] = c2[i][j];
      if ((i < minLevel) || (j < minLevel)) c1[i][j] = c2[i][j];
      c1[i][j] *= max_radius/scaling_factor;		/* Gebiet relativieren */
      }
   }

/* Initialisierung: belege den Block 
   (minLevel-1:minLevel,minLevel-1:minLevel) mit 1 */
nf = p*(1<<2*M);
init_pattern(T,nf,nf,20);
compute_bounding_boxes_gamma(&B1,W,E,p,M);
compute_bounding_boxes_square(&B2,W,E,p,M);
for (i=0; W[i].level<minLevel; i++)
{  for (j=0; j<=i; j++) set_pattern(T,i,j);
   }

/* explizite Rekursion */
for (ind1=0; W[ind1].level<M; ind1++)
{  m1 = W[ind1].level+1;
   
   /* bestimme aus den Bloecken (m1-1,1:m1-1) die Bloecke (m1,1:m1) */
   for (j=0; j<T->row_number[ind1]; j++)
   {  ind2 = T->index[ind1][j];
      m2 = W[ind2].level;
      if (m1 == m2+1)
      {  for (k=0; k<4; k++)
         {  if (wavelet_wavelet_criterion(B1,B2,W,E,W[ind1].son[k],ind2,c1[m1][m2],c2[m1][m2]))
            {  set_pattern(T,W[ind1].son[k],ind2);
	       for (l=0; l<4; l++)
	       {  if (wavelet_wavelet_criterion(B1,B2,W,E,W[ind1].son[k],W[ind2].son[l],c1[m1][m1],c2[m1][m1]))
	          {  set_pattern(T,W[ind1].son[k],W[ind2].son[l]);
	             }
	          }
	       }
	    if (wavelet_wavelet_criterion(B1,B2,W,E,W[ind2].son[k],ind1,c1[m1][m2],c2[m1][m2]))
            {  set_pattern(T,W[ind2].son[k],ind1);
	       for (l=0; l<4; l++)
	       {  if (wavelet_wavelet_criterion(B1,B2,W,E,W[ind2].son[k],W[ind1].son[l],c1[m1][m1],c2[m1][m1]))
	          {  set_pattern(T,W[ind2].son[k],W[ind1].son[l]);
	             }
	          }
	       }
	    }
	 }   
      else
      {  for (k=0; k<4; k++)
         {  if (wavelet_wavelet_criterion(B1,B2,W,E,W[ind1].son[k],ind2,c1[m1][m2],c2[m1][m2]))
            {  set_pattern(T,W[ind1].son[k],ind2);
	       }
	    }   
	 }   
      }
   }
   
/* Symmetrisierung der Steifigkeitsmatrix */
nnz = 0;
rn = (unsigned int*) calloc(nf,sizeof(unsigned int));
for (i=0; i<nf; i++)
{  rn[i] += T->row_number[i];
   for (j=0; j+1<T->row_number[i]; j++) rn[T->index[i][j]]++;
   }
for (i=0; i<nf; i++) 
{  T->index[i] = (unsigned int*) realloc(T->index[i],(rn[i]+1)*sizeof(unsigned int));
   T->max_row_number[i] = rn[i];
   T->index[i][rn[i]] = nf; 	/* setze Waechterelement */
   nnz += rn[i];
   }
for (i=0; i<nf; i++)
{  for (j=0; j+1<T->row_number[i]; j++) T->index[T->index[i][j]][T->row_number[T->index[i][j]]++] = i;
   }

/* Speicherplatz wieder freigeben */
/*printf("A-priori compression:            %.5f %%\n",100.0*nnz/nf/nf);*/
finish_pattern(T);
for (i=0; i<=M; i++)
{  free(c1[i]);
   free(c2[i]);
   }
free(rn);
free(c1);
free(c2);
free(B1);
free(B2);
return 100.0*nnz/nf/nf;
}
