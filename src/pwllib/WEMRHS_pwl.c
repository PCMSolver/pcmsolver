/**************
 *  WEMRHS_pwl.c  *
 **************/


/*=========================================================================*
 *  Bestimmt die rechte Seite im Wavelet-Galerkin-Verfahren fuer           *
 *  stueckweise lineare Wavelets bzgl. des modifizierten Skalarproduktes.  *
 *=========================================================================*/


#include <math.h>
#include <stdlib.h>
#include "constants.h"
#include "intvector_pwl.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "cubature.h"
#include "gauss_square.h"
#include "interpolate_pwl.h"
#include "WEMRHS_pwl.h"
#include "data.h"
#include "phi.h"


void WEMRHS_pwl1(rhs,W,E,T,p,M,nw)
/* testet die Neumann-Daten des gegebenen Potentials */
double		**rhs;		/* zu berechnende rechte Seite                */
wavelet		*W;		/* Waveletliste                               */
element		*E;		/* hierarchische Elementliste                 */
vector3		****T;		/* Oberflaecheninterpolation                  */
unsigned int	p;		/* Zahl der Patches                           */
unsigned int	M;		/* Zahl der Level                             */
unsigned int	nw;		/* Laenge von W				      */
{
unsigned int	N = 1 << M;	/* N*N Elemente pro Patch auf dem Level M     */
unsigned int	ne;		/* Anzahl der Elemente                        */
signed int	i, j;		/* Laufindizes durch die Wavelet/Elementliste */
double		h;		/* Schrittweite 			      */
cubature        *Q;		/* Kubatur-Formeln                            */
unsigned int	g = 0;		/* benoetigter Quadraturgrad                  */
double		c0, c1, c2, c3;	/* Werte der Integrale			      */
vector2		t;		/* Auswertepunkt der Gauss-Quadratur	      */
vector3		n_t;		/* Auswertung der Normale an den obigen Punkt */ 
unsigned int	k;	 	/* Laufindex fuer die Quadratur	              */
double		w;		/* Quadraturgewicht			      */
double		(*y)[4];	/* Array mit den Integralen (f,phi_k)         */
double		m;		/* Interpolationswert im Mittelpunkt	      */
double		e[4];		/* Interpolationswerte im Kantenmittelpunkt   */

/* Initialisierung */
ne = p*(4*N*N-1)/3;			/* Anzahl der Elemente */
init_Gauss_Square(&Q,g+1);		/* Kubatur-Formeln     */
y = (double (*)[4]) malloc(ne*sizeof(double[4]));
(*rhs) = (double*) malloc(nw*sizeof(double));

/* 1. Quadratur auf dem feinsten Level */
h = 1./N;
for (i=p*(N*N-1)/3; i<ne; i++)
{  c0 = c1 = c2 = c3 = 0;
   for (k=0; k<Q[g].nop; k++)
   {  t.x = h * (E[i].index_s + Q[g].xi[k].x);
      t.y = h * (E[i].index_t + Q[g].xi[k].y);
      n_t = n_Chi_pwl(t,T[E[i].patch],M);
      w = Q[g].w[k] * vector3_skalp(df(Chi_pwl(t,T[E[i].patch],M)),n_t);
      c0 += w * Phi0(Q[g].xi[k]);
      c1 += w * Phi1(Q[g].xi[k]);
      c2 += w * Phi2(Q[g].xi[k]);
      c3 += w * Phi3(Q[g].xi[k]);
      }
   y[i][0] = h * c0;
   y[i][1] = h * c1;
   y[i][2] = h * c2;
   y[i][3] = h * c3;
   }

/* 2. berechne Integrale der groeberen Level aus denen der feineren */
for (i=p*(N*N-1)/3-1; i>=0; i--)
{  m = 0.25 * (y[E[i].son[0]][2] + y[E[i].son[1]][3] + y[E[i].son[2]][0] + y[E[i].son[3]][1]);
   e[0] = 0.5 * (y[E[i].son[0]][1] + y[E[i].son[1]][0]);
   e[1] = 0.5 * (y[E[i].son[1]][2] + y[E[i].son[2]][1]);
   e[2] = 0.5 * (y[E[i].son[2]][3] + y[E[i].son[3]][2]);
   e[3] = 0.5 * (y[E[i].son[3]][0] + y[E[i].son[0]][3]);

   y[i][0] = 0.5 * (y[E[i].son[0]][0] + e[3] + e[0] + m);
   y[i][1] = 0.5 * (y[E[i].son[1]][1] + e[0] + e[1] + m);
   y[i][2] = 0.5 * (y[E[i].son[2]][2] + e[1] + e[2] + m);
   y[i][3] = 0.5 * (y[E[i].son[3]][3] + e[2] + e[3] + m);
   }

/* 3. setze Integrale (f,psi) zusammen */
for (i=0; i<nw; i++)
{  w = 0;
   for (j=0; j<W[i].element_number; j++)
   {  w += y[W[i].element[j]][0] * W[i].weight[j][0] \
         + y[W[i].element[j]][1] * W[i].weight[j][1] \
	 + y[W[i].element[j]][2] * W[i].weight[j][2] \
	 + y[W[i].element[j]][3] * W[i].weight[j][3];
      }
   (*rhs)[i] = w;
   }

/* Speicherplatz wieder freigeben */
free_Gauss_Square(&Q,g+1);
free(y);
return;
}


void WEMRHS_pwl2(rhs,W,E,T,p,M,nw)
/* testet die Dirichlet-Daten des gegebenen Potentials */
double		**rhs;		/* zu berechnende rechte Seite                */
wavelet		*W;		/* Waveletliste                               */
element		*E;		/* hierarchische Elementliste                 */
vector3		****T;		/* Oberflaecheninterpolation                  */
unsigned int	p;		/* Zahl der Patches                           */
unsigned int	M;		/* Zahl der Level                             */
unsigned int	nw;		/* Laenge von W				      */
{
unsigned int	N = 1 << M;	/* N*N Elemente pro Patch auf dem Level M     */
unsigned int	ne;		/* Anzahl der Elemente                        */
signed int	i, j;		/* Laufindizes durch die Wavelet/Elementliste */
double		h;		/* Schrittweite 			      */
cubature        *Q;		/* Kubatur-Formeln                            */
unsigned int	g = 0;		/* benoetigter Quadraturgrad                  */
double		c0, c1, c2, c3;	/* Werte der Integrale			      */
vector2		t;		/* Stuetzpunkt der Gauss-Quadratur auf Q      */
unsigned int	k;	 	/* Laufindex fuer die Quadratur	              */
double		w;		/* Quadraturgewicht			      */
double		(*y)[4];	/* Array mit den Integralen (f,phi_k)         */
double		m;		/* Interpolationswert im Mittelpunkt	      */
double		e[4];		/* Interpolationswerte im Kantenmittelpunkt   */

/* Initialisierung */
ne = p*(4*N*N-1)/3;			/* Anzahl der Elemente */
init_Gauss_Square(&Q,g+1);		/* Kubatur-Formeln     */
y = (double (*)[4]) malloc(ne*sizeof(double[4]));
(*rhs) = (double*) malloc(nw*sizeof(double));

/* 1. Quadratur auf dem feinsten Level */
h = 1./N;
for (i=p*(N*N-1)/3; i<ne; i++)
{  t.x = h * E[i].index_s;
   t.y = h * E[i].index_t;   
   c0 = 0.25*f(Chi_pwl(vector2_make(t.x        ,t.y        ),T[E[i].patch],M)); 
   c1 = 0.25*f(Chi_pwl(vector2_make(t.x+h-1e-14,t.y        ),T[E[i].patch],M));
   c2 = 0.25*f(Chi_pwl(vector2_make(t.x+h-1e-14,t.y+h-1e-14),T[E[i].patch],M));
   c3 = 0.25*f(Chi_pwl(vector2_make(t.x        ,t.y+h-1e-14),T[E[i].patch],M));
   c0 = c1 = c2 = c3 = 0;
   for (k=0; k<Q[g].nop; k++)
   {  t.x = h * (E[i].index_s + Q[g].xi[k].x);
      t.y = h * (E[i].index_t + Q[g].xi[k].y);
      w = Q[g].w[k] * f(Chi_pwl(t,T[E[i].patch],M));
      c0 += w * Phi0(Q[g].xi[k]);
      c1 += w * Phi1(Q[g].xi[k]);
      c2 += w * Phi2(Q[g].xi[k]);
      c3 += w * Phi3(Q[g].xi[k]);
      }
   y[i][0] = h * c0;
   y[i][1] = h * c1;
   y[i][2] = h * c2;
   y[i][3] = h * c3;
   }

/* 2. berechne Integrale der groeberen Level aus denen der feineren */
for (i=p*(N*N-1)/3-1; i>=0; i--)
{  m = 0.25 * (y[E[i].son[0]][2] + y[E[i].son[1]][3] + y[E[i].son[2]][0] + y[E[i].son[3]][1]);
   e[0] = 0.5 * (y[E[i].son[0]][1] + y[E[i].son[1]][0]);
   e[1] = 0.5 * (y[E[i].son[1]][2] + y[E[i].son[2]][1]);
   e[2] = 0.5 * (y[E[i].son[2]][3] + y[E[i].son[3]][2]);
   e[3] = 0.5 * (y[E[i].son[3]][0] + y[E[i].son[0]][3]);

   y[i][0] = 0.5 * (y[E[i].son[0]][0] + e[3] + e[0] + m);
   y[i][1] = 0.5 * (y[E[i].son[1]][1] + e[0] + e[1] + m);
   y[i][2] = 0.5 * (y[E[i].son[2]][2] + e[1] + e[2] + m);
   y[i][3] = 0.5 * (y[E[i].son[3]][3] + e[2] + e[3] + m);
   }

/* 3. setze Integrale (f,psi) zusammen */
for (i=0; i<nw; i++)
{  w = 0;
   for (j=0; j<W[i].element_number; j++)
   {  w += y[W[i].element[j]][0] * W[i].weight[j][0] \
         + y[W[i].element[j]][1] * W[i].weight[j][1] \
	 + y[W[i].element[j]][2] * W[i].weight[j][2] \
	 + y[W[i].element[j]][3] * W[i].weight[j][3];
      }
   (*rhs)[i] = w;
   }

/* Speicherplatz wieder freigeben */
free_Gauss_Square(&Q,g+1);
free(y);
return;
}

void WEMRHS_pwl2M(rhs,W,E,T,p,M,nw,potential,g)
/* testet die Dirichlet-Daten des gegebenen Potentials */
double		**rhs;		/* zu berechnende rechte Seite                */
wavelet		*W;		/* Waveletliste                               */
element		*E;		/* hierarchische Elementliste                 */
vector3		****T;		/* Oberflaecheninterpolation                  */
unsigned int	p;		/* Zahl der Patches                           */
unsigned int	M;		/* Zahl der Level                             */
unsigned int	nw;		/* Laenge von W				      */
unsigned int	g;		/* benoetigter Quadraturgrad                  */
double          *potential;
{
unsigned int	N = 1 << M;	/* N*N Elemente pro Patch auf dem Level M     */
unsigned int	ne;		/* Anzahl der Elemente                        */
signed int	i, j;		/* Laufindizes durch die Wavelet/Elementliste */
double		h;		/* Schrittweite 			      */
cubature        *Q;		/* Kubatur-Formeln                            */
double		c0, c1, c2, c3;	/* Werte der Integrale			      */
vector2		t;		/* Stuetzpunkt der Gauss-Quadratur auf Q      */
unsigned int	k;	 	/* Laufindex fuer die Quadratur	              */
double		w;		/* Quadraturgewicht			      */
double		(*y)[4];	/* Array mit den Integralen (f,phi_k)         */
double		m;		/* Interpolationswert im Mittelpunkt	      */
double		e[4];		/* Interpolationswerte im Kantenmittelpunkt   */
unsigned int    index;

/* Initialisierung */
ne = p*(4*N*N-1)/3;			/* Anzahl der Elemente */
init_Gauss_Square(&Q,g+1);		/* Kubatur-Formeln     */
y = (double (*)[4]) malloc(ne*sizeof(double[4]));
(*rhs) = (double*) malloc(nw*sizeof(double));

/* 1. Quadratur auf dem feinsten Level */
h = 1./N;
for (i=p*(N*N-1)/3; i<ne; i++)
{  t.x = h * E[i].index_s;
   t.y = h * E[i].index_t;   
   c0 = 0.25*f(Chi_pwl(vector2_make(t.x        ,t.y        ),T[E[i].patch],M)); 
   c1 = 0.25*f(Chi_pwl(vector2_make(t.x+h-1e-14,t.y        ),T[E[i].patch],M));
   c2 = 0.25*f(Chi_pwl(vector2_make(t.x+h-1e-14,t.y+h-1e-14),T[E[i].patch],M));
   c3 = 0.25*f(Chi_pwl(vector2_make(t.x        ,t.y+h-1e-14),T[E[i].patch],M));
   c0 = c1 = c2 = c3 = 0;
   for (k=0; k<Q[g].nop; k++)
   { 
     /*
     t.x = h * (E[i].index_s + Q[g].xi[k].x);
     t.y = h * (E[i].index_t + Q[g].xi[k].y);
     w = Q[g].w[k] * f(Chi_pwl(t,T[E[i].patch],M));
     */
     index = (E[i].patch*N*N + E[i].index_t*N + E[i].index_s ) * Q[g].nop + k;
     w = Q[g].w[k] * potential[index];
     c0 += w * Phi0(Q[g].xi[k]);
     c1 += w * Phi1(Q[g].xi[k]);
     c2 += w * Phi2(Q[g].xi[k]);
     c3 += w * Phi3(Q[g].xi[k]);
      }
   y[i][0] = h * c0;
   y[i][1] = h * c1;
   y[i][2] = h * c2;
   y[i][3] = h * c3;
   }

/* 2. berechne Integrale der groeberen Level aus denen der feineren */
for (i=p*(N*N-1)/3-1; i>=0; i--)
{  m = 0.25 * (y[E[i].son[0]][2] + y[E[i].son[1]][3] + y[E[i].son[2]][0] + y[E[i].son[3]][1]);
   e[0] = 0.5 * (y[E[i].son[0]][1] + y[E[i].son[1]][0]);
   e[1] = 0.5 * (y[E[i].son[1]][2] + y[E[i].son[2]][1]);
   e[2] = 0.5 * (y[E[i].son[2]][3] + y[E[i].son[3]][2]);
   e[3] = 0.5 * (y[E[i].son[3]][0] + y[E[i].son[0]][3]);

   y[i][0] = 0.5 * (y[E[i].son[0]][0] + e[3] + e[0] + m);
   y[i][1] = 0.5 * (y[E[i].son[1]][1] + e[0] + e[1] + m);
   y[i][2] = 0.5 * (y[E[i].son[2]][2] + e[1] + e[2] + m);
   y[i][3] = 0.5 * (y[E[i].son[3]][3] + e[2] + e[3] + m);
   }

/* 3. setze Integrale (f,psi) zusammen */
for (i=0; i<nw; i++)
{  w = 0;
   for (j=0; j<W[i].element_number; j++)
   {  w += y[W[i].element[j]][0] * W[i].weight[j][0] \
         + y[W[i].element[j]][1] * W[i].weight[j][1] \
	 + y[W[i].element[j]][2] * W[i].weight[j][2] \
	 + y[W[i].element[j]][3] * W[i].weight[j][3];
      }
   (*rhs)[i] = w;
   }

/* Speicherplatz wieder freigeben */
free_Gauss_Square(&Q,g+1);
free(y);
return;
}
