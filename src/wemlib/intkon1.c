/***************
 *  IntKon1.c  *
 ***************/


/*===========================================================*
 *  Fernfeld-Quadratur-Routine:			             *
 *  Die in der Funktion init_randwerte vorab berechneten     *
 *  Auswertepunkte und Gewichte der Gauss-Quadratur werden   *	
 *  in IntKon1 zum entsprechenden Integral zusammengefuegt.  *
 *===========================================================*/
 

#include <stdlib.h>
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "kern.h"
#include "cubature.h"
#include "gauss_square.h"
#include "interpolate.h"
#include "intkon1.h"
#include "integrate.h"


void init_randwerte(RW,g_max)
/* Initialisiert die Randwerte */
randwerte	**RW;
unsigned int	g_max;
{
unsigned int	i;

(*RW) = (randwerte*) malloc(g_max*sizeof(randwerte));
for (i=0; i<g_max; i++) 
{  (*RW)[i].Chi   = (vector3*) malloc((i+1)*(i+1)*sizeof(vector3));
   (*RW)[i].n_Chi = (vector3*) malloc((i+1)*(i+1)*sizeof(vector3));
   (*RW)[i].det_dChi = (double*) malloc((i+1)*(i+1)*sizeof(double));
   }
return;
}


void reset_randwerte(RW,g_max)
/* Resetted die Randwerte */
randwerte	*RW;
unsigned int	g_max;
{
unsigned int	i;
for (i=0; i<g_max; i++) RW[i].nop = 0;
return;
}


void free_randwerte(RW,g_max)
/* Gibt den Speicherplatz der Randwerte frei */
randwerte	**RW;
unsigned int	g_max;
{
unsigned int	i;

for (i=0; i<g_max; i++) 
{  free((*RW)[i].Chi);
   free((*RW)[i].n_Chi);
   free((*RW)[i].det_dChi);
   }
free(*RW);
return;
}


void IntKon1(c,element1,element2,RW,Q1,Q2,P,M,SingleLayer,DoubleLayer)
/* No-Problem-Quadrature-Routine -> modifiziertes Skalarprodukt */
double		*c;
element		*element1, *element2;
randwerte	*RW;
cubature	*Q1, *Q2;
vector3		****P;
unsigned int	M;
double		SingleLayer(), DoubleLayer();
{
double		h_s;
vector2		xi, eta;
vector3		y, n_y;
unsigned int	i, j;
double		h_t = 1./(1 << element2->level);

/* falls noetig, berechne Randwerte */
if (RW->nop == 0) 
{  RW->nop = Q1->nop;
   h_s = 1./(1 << element1->level);
   for (i=0; i<Q1->nop; i++)
   {  xi.x = h_s*(element1->index_s+Q1->xi[i].x);
      xi.y = h_s*(element1->index_t+Q1->xi[i].y);
      RW->Chi[i] = Chi(xi,P[element1->patch],M);
      RW->n_Chi[i] = n_Chi(xi,P[element1->patch],M);
      RW->det_dChi[i] = h_s * Q1->w[i];
      }
   }

/* Quadratur */
c[0] = c[1] = c[2] = 0;
for (i=0; i<Q2->nop; i++)
{  eta.x = h_t * (element2->index_s + Q2->xi[i].x);
   eta.y = h_t * (element2->index_t + Q2->xi[i].y);
   y = Chi(eta,P[element2->patch],M);
   n_y = n_Chi(eta,P[element2->patch],M);
   for (j=0; j<RW->nop; j++) 
   {  c[0] += Q2->w[i] * RW->det_dChi[j] * SingleLayer(RW->Chi[j],y);
      c[1] += Q2->w[i] * RW->det_dChi[j] * DoubleLayer(RW->Chi[j],y,n_y);
      c[2] += Q2->w[i] * RW->det_dChi[j] * DoubleLayer(y,RW->Chi[j],RW->n_Chi[j]);
      }
   }

c[0] *= h_t;		/* L^2-normiert */
c[1] *= h_t;
c[2] *= h_t;
return;
}
