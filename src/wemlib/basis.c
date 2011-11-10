/*************
 *  Basis.c  *
 *************/
 

/*===================================================*
 *  Erzeugt die fuer das Wavelet-Galerkin-Verfahren  *
 *  notwendige Liste der Elemente ueber alle Gitter  *
 *  und die Waveletliste.                            *
 *===================================================*/
 
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "sparse.h"
#include "basis.h"
#include "mask.h"


/*=================================================================*
 *  D E K L A R A T I O N   D E R   H I L F S F U N K T I O N E N  *
 *=================================================================*/

void add_element(wavelet *w, element *E, double weight, unsigned int index);
/* fuegt zum Wavelet w das mit weight gewichtete Element index hinzu */


/*===============================================================*
 *  D E F I N I T I O N   D E R   H I L F S F U N K T I O N E N  *
 *===============================================================*/

void unify(d,r,d1,r1,d2,r2)
/* bildet die Vereinigung K(d,r) = K(d1,r1) \cup K(d2,r2) */
vector3		*d, d1, d2;
double		*r, r1, r2;
{
vector3		z;
double		norm;

z.x = d1.x-d2.x;
z.y = d1.y-d2.y;
z.z = d1.z-d2.z;
norm = sqrt(z.x*z.x+z.y*z.y+z.z*z.z);

if (norm+r2 <= r1) 		/* K(d2,r2) \subset K(d1,r1) */
{  *d = d1;
   *r = r1;
   }
else if (norm+r1 <= r2)		/* K(d1,r1) \subset K(d2,r2) */
{  *d = d2;
   *r = r2;
   }
else		/* die Vereinigungsmenge ist keine Kugel */
{  (*d).x = 0.5*(d1.x+d2.x+(r1-r2)/norm*z.x);
   (*d).y = 0.5*(d1.y+d2.y+(r1-r2)/norm*z.y);
   (*d).z = 0.5*(d1.z+d2.z+(r1-r2)/norm*z.z);
   *r = 0.5*(r1+r2+norm);
   }
return;
}


void add_element(w,E,weight,index)
/* fuegt zum Wavelet w das mit weight gewichtete Element index hinzu */
wavelet		*w;		/* gegebenes Wavelet                              */
element		*E;		/* hierarchische Elementliste                     */
double		weight;		/* Gewicht der zu addierenden Skalierungsfunktion */
unsigned int	index;		/* Index des zu addierenden Elements              */
{
if (weight == 0) return;	/* es ist nix zu tun */

if (w->element_number%delta == 0)
{  w->element = (unsigned int*) realloc(w->element,(w->element_number+delta)*sizeof(unsigned int));
   w->weight  = (double*)       realloc(w->weight, (w->element_number+delta)*sizeof(double));
   }
   
w->element[w->element_number] = index;
w->weight[w->element_number] = weight;
w->element_number++;
return;
}


/*===========================*
 *  E L E M E N T L I S T E  *
 *===========================*/

unsigned int generate_elementlist(E,P,F,p,M)
/* erstellt die hierarchsische Elementliste E */
element		**E;		/* hierarchische Elementliste             */
vector3		*P;		/* Punkteliste der Einskalenbasis         */
unsigned int	**F;		/* Elementliste der Einskalenbasis        */
unsigned int	p;		/* Anzahl der Patches                     */
unsigned int	M;		/* 2^M*2^M Elemente pro Patch             */
{
signed int	m;		/* Laufindex fuer das Level               */
unsigned int	n;		/* p*n*n Elemente auf Level m             */
unsigned int	N = 1 << M;	/* p*N*N Elemente auf Level M             */
unsigned int    i1, i2, i3;	/* Laufindizes durch Gitter m             */
unsigned int	zf, ze;		/* Zaehler fuer Elementlisten             */
unsigned int	ne;		/* Laenge von E                           */
vector3		d1, d2;		/* Umkreismittelpunkt von Element 0 und 2 */
double		r1, r2;		/* Umkreisradius      von Element 1 und 3 */

ne = p*(4*N*N-1)/3;	/* Laenge von E */
(*E) = (element*) calloc(ne,sizeof(element));

/* initialisiere feinstes Gitter */
zf = 0;
ze = p*(N*N-1)/3;
for (i1=0; i1<p; i1++)
{  for (i2=0; i2<N; i2++)
   {  for (i3=0; i3<N; i3++)
      {  
         (*E)[ze].level = M;			/* Level des Elements        */
	 (*E)[ze].patch = i1;			/* Patch des Elements        */
	 (*E)[ze].index_s = i3;			/* Index in s-Richtung       */
	 (*E)[ze].index_t = i2;			/* Index in t-Richtung       */
	 (*E)[ze].vertex[0] = F[zf][0];		/* 1. Eckpunkt des Elements  */
	 (*E)[ze].vertex[1] = F[zf][1];		/* 2. Eckpunkt des Elements  */
	 (*E)[ze].vertex[2] = F[zf][2];		/* 3. Eckpunkt des Elements  */
	 (*E)[ze].vertex[3] = F[zf][3];		/* 4. Eckpunkt des Elements  */

	 /* setze die Soehne des Vaterelements */
	 (*E)[ze].son[0] = (*E)[ze].son[1] = (*E)[ze].son[2] = (*E)[ze].son[3] = ne;

	 /* bestimme Mittelpunkt und Radius des Elementumkreises */
	 unify(&d1,&r1,P[F[zf][0]],0,P[F[zf][2]],0);
	 unify(&d2,&r2,P[F[zf][1]],0,P[F[zf][3]],0);
	 unify(&(*E)[ze].midpoint,&(*E)[ze].radius,d1,r1,d2,r2);

	 ze++;
	 zf++;
	 }
      }
   }

/* Schleife ueber die groeberen Gitter */
ze = ne;
for (m=M-1; m>=0; m--)
{  n = 1 << m;
   ze -= 5*p*n*n;
   zf = ze+p*n*n;
   for (i1=0; i1<p; i1++)
   {  for (i2=0; i2<n; i2++)
      {  for (i3=0; i3<n; i3++)
         {  
	    (*E)[ze].level = m;					/* Level des Elements       */
	    (*E)[ze].patch = i1;				/* Patch des Elements       */
	    (*E)[ze].index_s = i3;				/* Index in s-Richtung      */
	    (*E)[ze].index_t = i2;				/* Index in t-Richtung      */
	    (*E)[ze].vertex[0] = (*E)[zf      ].vertex[0];	/* 1. Eckpunkt des Elements */
	    (*E)[ze].vertex[1] = (*E)[zf    +1].vertex[1];	/* 2. Eckpunkt des Elements */
	    (*E)[ze].vertex[2] = (*E)[zf+2*n+1].vertex[2];	/* 3. Eckpunkt des Elements */
	    (*E)[ze].vertex[3] = (*E)[zf+2*n  ].vertex[3];	/* 4. Eckpunkt des Elements */
	    (*E)[ze].son[0] = zf;				/* 1. Sohn des Elements     */
	    (*E)[ze].son[1] = zf+1;				/* 2. Sohn des Elements     */
	    (*E)[ze].son[2] = zf+2*n+1;				/* 3. Sohn des Elements     */
	    (*E)[ze].son[3] = zf+2*n;				/* 4. Sohn des Elements     */

      	    /* setze in den Soehnen das Vaterelement */
	    (*E)[zf].father = (*E)[zf+1].father = (*E)[zf+2*n].father = (*E)[zf+2*n+1].father = ze;
	    
	    /* bestimme Mittelpunkt und Radius des Elementumkreises */
	    unify(&d1,&r1,(*E)[zf  ].midpoint,(*E)[zf  ].radius,(*E)[zf+2*n+1].midpoint,(*E)[zf+2*n+1].radius);
	    unify(&d2,&r2,(*E)[zf+1].midpoint,(*E)[zf+1].radius,(*E)[zf+2*n  ].midpoint,(*E)[zf+2*n  ].radius);
	    unify(&(*E)[ze].midpoint,&(*E)[ze].radius,d1,r1,d2,r2);

	    zf += 2;
	    ze++;
	    }
	 zf += 2*n;
	 }
      }
   }
return(ne);
}


void complete_elementlist(W,E,p,M)
/* erstellt die Liste zum Zugriff auf die Wavelets */
wavelet		*W;		/* Liste der Wavelets                 */
element		*E;		/* hierarchische Elementliste         */
unsigned int	p;		/* Anzahl der Patches                 */
unsigned int	M;		/* 2^M*2^M Elemente pro Patch         */
{
unsigned int	N = 1 << M;	/* p*N*N Elemente auf Level M         */
unsigned int	nw;		/* Laenge der Waveletliste            */
unsigned int	ne;		/* Laenge hierarchische Elementliste  */
unsigned int    i, j;		/* Laufindizes durch die Waveletliste */
unsigned int    k;		/* zu bearbeitendes Element           */

nw = p*N*N;		/* Laenge von W */
ne = p*(4*N*N-1)/3;;	/* Laenge von E */

/* Elementliste updaten */
for (i=0; i<nw; i++)
{  for (j=0; j<W[i].element_number; j++)
   {  k = W[i].element[j];
      if (E[k].wavelet_number%delta == 0)
      {  E[k].wavelet = (unsigned int*) realloc(E[k].wavelet,(E[k].wavelet_number+delta)*sizeof(unsigned int));
         }
      E[k].wavelet[E[k].wavelet_number++] = i;
      }
   }

/* Speicherplatz der Elementliste optimieren */
for (i=0; i<ne; i++) E[i].wavelet = (unsigned int*) realloc(E[i].wavelet,E[i].wavelet_number*sizeof(unsigned int));
return;
}


void free_elementlist(E,p,M)
/* gibt den Speicherplatz der hierarchischen Elementliste E frei */
element		**E;		/* hierarchische Elementliste */
unsigned int	p;		/* Anzahl der Patches         */
unsigned int	M;		/* 2^M*2^M Elemente pro Patch */
{
unsigned int	ne;		/* Laenge von E               */
unsigned int	i;		/* Laufindex durch E          */

ne = p*(4*(1<<2*M)-1)/3;	/* Anzahl aller Elemente */

for (i=0; i<ne; i++) free((*E)[i].wavelet);
free(*E);
}


/*===========================*
 *  W A V E L E T L I S T E  *
 *===========================*/

unsigned int generate_waveletlist(W,E,p,M)
/* erstellt die Waveletliste W */
wavelet		**W;		/* Liste der Wavelets                     */
element		*E;		/* hierarchische Elementliste             */
unsigned int	p;		/* Anzahl der Patches                     */
unsigned int	M;		/* 2^M*2^M Elemente pro Patch             */
{
unsigned int	nw;		/* Anzahl der Wavelets                    */
unsigned int	N = 1 << M;	/* p*N*N Elemente auf Level M             */
unsigned int	n;		/* p*n*n Elemente auf Level m             */
signed int	m;		/* Laufindex fuer das Level               */
unsigned int    i1, i2, i3;	/* Laufindizes durch Gitter m             */
unsigned int	zw;		/* Nullwavelet von Patch i1 auf Level m   */
unsigned int	ze;		/* Nullelement von Patch i1 auf Level m   */
sparse		T, L;		/* Maskenmatrizen                         */
unsigned int	s, t;		/* Laufindizes durch die Maskenmatrizen   */
wavelet		*w;		/* Zeiger auf das zu bearbeitende Wavelet */

/* SCHLEIFE UEBER DIE GITTER */
nw = p*N*N;		/* Laenge der Waveletliste */
(*W) = (wavelet*) calloc(nw,sizeof(wavelet));
for (m=M; m>=minLevel; m--)
{
   dwt_mask(&T,&L,m);	/* berechne die Masken T und L  */
   n = 1 << (m-1);	/* p*n*n Elemente auf Level m-1 */

   /* 1. bilde Skalierungsfunktionen und Grund-Wavelets */
   for (i1=0; i1<p; i1++)
   {  zw = i1*n*n;		/* Nullpunkt von Patch i1 auf Level m-1 in der Waveletliste */
      ze = p*(4*n*n-1)/3+4*zw;	/* Nullpunkt von Patch i1 auf Level m   in der Elementliste */

      for (i2=0; i2<n; i2++)
      {  for (i3=0; i3<n; i3++)
         {  
	    /* 1. psi(s) * phi(t) */
	    w = &(*W)[zw+p*n*n+n*i2+i3];
	    w->level = m;				/* Level des Wavelets */
	    w->son[0] = 4*p*n*n+4*zw+4*i2*n+2*i3;	/* bestimme Soehne    */
	    w->son[1] = w->son[0] + 1;
	    w->son[2] = w->son[0] + 2*n;
	    w->son[3] = w->son[0] + 2*n + 1;
	    for (s=0; s<T.row_number[i3]; s++)	/* addiere Skalierungsfunktionen von Level m */
	    {  add_element(w,E,0.5*T.value[i3][s],ze+(4*n*i2    )+T.index[i3][s]);
	       add_element(w,E,0.5*T.value[i3][s],ze+(4*n*i2+2*n)+T.index[i3][s]);
               }

	    /* 2. phi_1(s) * psi(t) */
	    w = &(*W)[zw+2*p*n*n+n*i2+i3];
	    w->level = m;				/* Level des Wavelets */
	    w->son[0] = 8*p*n*n+4*zw+4*i2*n+2*i3;	/* bestimme Soehne    */
	    w->son[1] = w->son[0] + 2*n;
	    w->son[2] = w->son[0] + 4*p*n*n;
	    w->son[3] = w->son[0] + 4*p*n*n + 2*n;
	    for (t=0; t<T.row_number[i2]; t++)
	    {  add_element(w,E,0.5*T.value[i2][t],ze+(2*n*T.index[i2][t])+(2*i3));
               }

	    /* 3. phi_2(s) * psi(t) */
	    w = &(*W)[zw+3*p*n*n+n*i2+i3];
	    w->level = m;				/* Level des Wavelets */
	    w->son[0] = 8*p*n*n+4*zw+4*i2*n+2*i3+1;	/* bestimme Soehne    */
	    w->son[1] = w->son[0] + 2*n;
	    w->son[2] = w->son[0] + 4*p*n*n;
	    w->son[3] = w->son[0] + 4*p*n*n + 2*n;
	    for (t=0; t<T.row_number[i2]; t++)
	    {  add_element(w,E,0.5*T.value[i2][t],ze+(2*n*T.index[i2][t])+(2*i3+1));
               }
       	    }
	 }
      }
   free_sparse(&T);
   }

/* noch Skalierungsfunktionen (Level minLevel-1) hinzufuegen */
zw = 0;
n = 1 << (minLevel-1);
for (i1=0; i1<p; i1++)
{  for (i2=0; i2<n; i2++)
   {  for (i3=0; i3<n; i3++)
      {  (*W)[zw].level = minLevel-1;		/* Level der Skalierungsfunktion  */
	 (*W)[zw].son[0] = zw;			/* Soehne der Skalierungsfunktion */
	 (*W)[zw].son[1] = zw+p*n*n;
	 (*W)[zw].son[2] = zw+2*p*n*n;
	 (*W)[zw].son[3] = zw+3*p*n*n;
	 add_element(&(*W)[zw],E,1,zw);		/* Element + Einschliessung */
	 zw++;
	 }
      }
   }
return(nw);
}


void set_quadrature_level(W,E,p,M)
/* verfeinert Grobgitterelemente */
wavelet		*W;			/* Liste der Wavelets                          */
element		*E;			/* hierarchische Elementliste                  */
unsigned int	p;			/* Anzahl der Patches                          */
unsigned int	M;			/* 2^M*2^M Elemente pro Patch                  */
{
unsigned int	nw;			/* Laenge von W                                */
unsigned int	i, j, k, l;		/* Laufindizes durch Wavelet/Elementliste      */
unsigned int	ind;			/* Index des zu untersuchenden Elements        */
unsigned int	minLevel;		/* Minimales Level fuer die Quadratur          */
unsigned int	noe;			/* number of elements eines Wavelets           */

nw = p*(1<<M)*(1<<M);			/* Anzahl der Wavelets            */
minLevel = min_quadrature_level;	/* setze minimales Quadraturlevel */
if (minLevel > M) minLevel = M;		/* sonst gibt's Aerger            */

for (i=0; (i<nw)&&(W[i].level<=minLevel); i++)
{  for (j=W[i].level; (j<=minLevel); j++)
   {  noe = W[i].element_number;
      for (k=0; k<noe; k++)
      {  ind = W[i].element[k];	           /* zu untersuchendes Element  */
         if (E[ind].level < minLevel)      /* also zu grob -> verfeinere */ 
         {  W[i].element[k] = E[ind].son[0];
            W[i].weight[k] *= 0.5;
	    W[i].element = (unsigned int*) realloc(W[i].element,(W[i].element_number+3)*sizeof(unsigned int));
	    W[i].weight  = (double*)       realloc(W[i].weight ,(W[i].element_number+3)*sizeof(double));
	    for (l=1; l<4; l++)
	    {  W[i].element[W[i].element_number] = E[ind].son[l];
	       W[i].weight[W[i].element_number]  = W[i].weight[k];
	       W[i].element_number++;
	       }
	    }
         }
      }
   }
return;
}


void simplify_waveletlist(W,E,p,M)
/* optimiert die Waveletliste W */
wavelet		*W;			/* Liste der Wavelets                          */
element		*E;			/* hierarchische Elementliste                  */
unsigned int	p;			/* Anzahl der Patches                          */
unsigned int	M;			/* 2^M*2^M Elemente pro Patch                  */
{
unsigned int	nw;			/* Laenge von W                                */
unsigned int	i, k, l;		/* Laufindizes durch Wavelet/Elementliste      */
signed int	j;			/* Laufindex durch Wavelet/Elementliste        */
unsigned int	ind;			/* Index des zu untersuchenden Elements        */
unsigned int	s0, s1, s2, s3;		/* Laufindizes fuer die Suche nach den Soehnen */
unsigned int	noe;			/* number of elements eines Wavelets           */
unsigned int	*prototype;		/* Liste fuer den Zugriff Prototyp -> Wavelet  */
unsigned int	prototype_number;	/* Anzahl der diversen Protptypen	       */
unsigned int	minLevel;		/* Minimales Level fuer die Quadratur          */

nw = p*(1<<M)*(1<<M);			/* Anzahl der Wavelets            */
minLevel = min_quadrature_level;	/* setze minimales Quadraturlevel */

/* 1. vereinfache die Wavelets */
for (i=0; i<nw; i++)
{  if (W[i].level > minLevel)
   {  /* Untersuche, ob Feingitterelemente durch ihren Vater ersetzt werden   */
      /* koennen. Falls ja, tue dies, und setze das Gewicht der Soehne auf 0. */
      noe = W[i].element_number;	/* Zahl der Eintraege in der Elementliste des Wavelets */
      for (s0=0; s0<noe; s0++)
      {  ind = W[i].element[s0];	   /* zu untersuchendes Element      		       */
         if (E[ind].level == W[i].level)   /* also Feingitterelement -> suche nach den Soehnen */
         {  /* ueberpruefe, ob zu untersuchendes Element der 0. Sohn ist */
	    ind = E[ind].father;		/* Vaterelement */
	    if (W[i].element[s0] == E[ind].son[0])
            {  /* 0. Sohn vorhanden -> ueberpruefe auf 1. Sohn */
	       for (s1=0; (s1<noe) && (E[ind].son[1]!=W[i].element[s1]); s1++);
	       if ((s1 < noe) && (W[i].weight[s1] == W[i].weight[s0]))
	       {  /* 0.+1. Sohn vorhanden -> ueberpruefe auf 2. Sohn */
	          for (s2=0; (s2<noe) && (E[ind].son[2]!=W[i].element[s2]); s2++);
	          if ((s2 < noe) && (W[i].weight[s2] == W[i].weight[s0]))
	          {  /* 0.+1.+2. Sohn vorhanden -> ueberpruefe auf 3. Sohn */
	             for (s3=0; (s3<noe) && (E[ind].son[3]!=W[i].element[s3]); s3++);
	             if ((s3 < noe) && (W[i].weight[s3] == W[i].weight[s0]))
	             {  /* alle 4 Soehne sind vorhanden und haben gleiche Gewichte */
	                W[i].element[s0] = ind;
		        W[i].weight[s0] *= 2;
                        W[i].weight[s1] = W[i].weight[s2] = W[i].weight[s3] = 0;
		        }
		     }
		  }
	       }
	    }
         }

      /* Groesse der Elementliste anpassen */
      k = 0;
      while (k < W[i].element_number)
      {  if (W[i].weight[k] == 0)
         {  W[i].element_number--;		/* streiche Element */
            for (l=k; l<W[i].element_number; l++)
	    {  W[i].element[l] = W[i].element[l+1];
	       W[i].weight[l] = W[i].weight[l+1];
	       }
	    }
         else k++;
         }      
      W[i].element = (unsigned int*) realloc(W[i].element,W[i].element_number*sizeof(unsigned int));
      W[i].weight  = (double*)       realloc(W[i].weight, W[i].element_number*sizeof(double));
      }
   }

/* 2. bestimme die Liste der Prototypen */
prototype = NULL;
prototype_number = 0;

for (i=0; i<nw; i++)
{  for (j=prototype_number-1; j>=0; j--)
   {
      if (W[i].element_number == W[prototype[j]].element_number)
      {  
      	 /* gleiche Elementanzahl: gleiche Gewichte? */
         for (k=0; (k<W[i].element_number) && (W[i].weight[k] == W[prototype[j]].weight[k]); k++);

	 if (k == W[i].element_number) 			/* alle Gewichte gleich ==> ersetze */
	 {  free(W[i].weight);	        		/* die Liste der Gewichte durch die */
	    W[i].weight = W[prototype[j]].weight;	/* des Prototypen und verlasse die  */
	    break;					/* Schleife bzgl. j.                */
	    }
	 }
      }

   if (j == -1)		/* keinen passenden Prototyp gefunden */
   {  if (prototype_number%delta == 0) prototype = (unsigned int*) realloc(prototype,(prototype_number+delta)*sizeof(unsigned int));
      prototype[prototype_number++] = i;
      }
   }

free(prototype);
/*printf("%d prototypes\n",prototype_number);*/
return;
}


void free_waveletlist(W,p,M)
/* gibt den Speicherplatz der Waveletliste W frei */
wavelet		**W;			/* Liste der Wavelets                         */
unsigned int	p;			/* Anzahl der Patches                         */
unsigned int	M;			/* 2^M*2^M Elemente pro Patch 		      */
{
unsigned int	nw;			/* Laenge von W               		      */
unsigned int	i;			/* Laufindex durch W          		      */
unsigned int	*prototype;		/* Liste fuer den Zugriff Prototyp -> Wavelet */
unsigned int	prototype_number;	/* Anzahl der diversen Protptypen	      */
signed int	j;			/* Laufindex durch die Liste der Prototypen   */

prototype = NULL;
prototype_number = 0;
nw = p*(1<<M)*(1<<M);	/* Anzahl der Wavelets */

for (i=0; i<nw; i++)
{  free((*W)[i].element);

   /* suche nach passendem Prototypen */
   for (j=prototype_number-1; j>=0; j--)
   {  if ((*W)[i].weight == (*W)[prototype[j]].weight) break;
      }

   if (j == -1)		/* keinen passenden Prototyp gefunden */
   {  if (prototype_number%delta == 0) prototype = (unsigned int*) realloc(prototype,(prototype_number+delta)*sizeof(unsigned int));
      prototype[prototype_number++] = i;
      }
   }

for (i=0; i<prototype_number; i++) free((*W)[prototype[i]].weight);
free(prototype);
free(*W);
}


/*===================================*
 *  A B S T A N D S F U N K T I O N  *
 *===================================*/
 
double distance(element1,element2)
/* Berechnet den Abstand zwischen den Elementen element1 und element2 */
element 	*element1, *element2;	/* Pointer auf die zwei Elemente */
{
double		dx, dy, dz;		/* x/y/z-Abstand zweier Elemente */
double		c1, c2;			/* Skalierungsfaktoren           */
double		dist;			/* zu berechnender Abstand       */

dx = element1->midpoint.x-element2->midpoint.x;
dy = element1->midpoint.y-element2->midpoint.y;
dz = element1->midpoint.z-element2->midpoint.z;
dist = sqrt(dx*dx+dy*dy+dz*dz) - element1->radius - element2->radius;
c1 = element1->radius*(1<<element1->level);
c2 = element2->radius*(1<<element2->level);

if (c1 < c2) return(scaling_factor*dist/c2); else return(scaling_factor*dist/c1);
}
