/*****************
 *  Integrate.c  *
 *****************/
 
 
/*=============================================================*
 *  Enthaelt alle Routinen, die im Wavelet-Galerkin-Verfahren  *
 *  fuer das aufsplitten der Integrale benoetigt werden.       *
 *=============================================================*/
 

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "intvector_pwl.h"
#include "vector2.h"
#include "vector3.h"
#include "basis_pwl.h"
#include "trafos_pwl.h"
#include "cubature.h"
#include "intlin1.h"
#include "intlin2.h"
#include "intlin3.h"
#include "intlin4.h"
#include "integrate_pwl.h"


void element_element_interaction_pwl(c,P,E,ind1,ind2,RW,Q,R,M,prec,SingleLayer,
                                     DoubleLayer,Identity)
/* Zerlegungsalgorithmus fuer die Integration Element ind1 mit Element ind2 */
double		*c;		/* zu berechnende Integrale                      */
vector3		*P;		/* Punkteliste                                   */
element_pwl		*E;		/* hierarchische Elemntliste                     */
unsigned int	ind1, ind2;	/* Indizes der Integranden                       */
randwerte	*RW;		/* Randwerte bezueglich Element ind1             */
cubature	*Q; 		/* Kubatur-Formeln                               */
vector3		****R;		/* Koeffizienten zur Oberflaecheninterpolation   */
unsigned int	M;		/* Zahl der Level                                */
double		prec;		/* Quadraturgenauigkeit                          */
double		SingleLayer();	/* kernel of the single layer operator           */ 
double		DoubleLayer();  /* kernel of the double layer operator           */
double          Identity;       /* correctly scaled identity operator            */
{
signed int	j;		/* Listenindex eines Elements der Integralmatrix */
signed int	g1, g2;		/* benoetigte Quadraturgrade                     */
double		dist;		/* Abstand zwischen den Elementen ind1 und ind2  */
unsigned int	s, t;		/* Indizes der Drehungen bei sing. Integralen    */
unsigned int	CASE;		/* Fallunterscheidung bei Integration            */
double		a[192];		/* Werte der Integrale bei Unterteilung          */
unsigned int	k;		/* Laufindex fuer die Interpolation              */
double		m;		/* Interpolationswert im Mittelpunkt		 */
double		e[4];		/* Interpolationswerte im Kantenmittelpunkt	 */

/* untersuche zuerst, ob das Integral schon berechnet worden ist */
if (ind1 >= ind2) 
{  j = search_integral_pwl(&E[ind1].interaction,ind2);
   if (j != -1) 
   {  memcpy(c,E[ind1].interaction.value[j].sub,48*sizeof(double));
      return;
      }
   }
else 
{  j = search_integral_pwl(&E[ind2].interaction,ind1);
   if (j != -1) 
   {  permutate(c,E[ind2].interaction.value[j].sub);
      return;
      }
   }

/* bestimme den Abstand zwischen den Elementen */
dist = distance_pwl(&E[ind1],&E[ind2]);

/* Quadratur mit Genauigkeit prec */
if (E[ind1].level == E[ind2].level)	/* nicht unterteilen, da beide Elemente auf gleichem Level */
{  
  /* Falleinteilung fuer Quadratur */
   if (dist > eps)         CASE = 1;  		  		/* keine Gemeinsamkeiten   */
   else if (ind1 == ind2)  CASE = 2;  				/* gleiches Patch          */ 
   else CASE = compare_pwl(P,E[ind1].vertex,E[ind2].vertex,&s,&t);	/* gemeinsame Kante/Punkt? */

   /* Quadratur mit Genauigkeit prec */
   if (dist*(1<<E[ind1].level) < 1) dist = 1./(1<<E[ind1].level);
   quadrature_grade_pwl(&g1,&g2,E[ind1].level,E[ind2].level,dist,prec);

   /* Wahl der Quadraturroutine aufgrund der Falleinteilung */
   switch (CASE)
   {  case 1:  IntLin1(c,&E[ind1],&E[ind2],&RW[g1],&Q[g1],&Q[g1],R,M,SingleLayer,DoubleLayer);
	       break;
      case 2:  IntLin2(c,&E[ind1],&Q[g1],R,M,SingleLayer,DoubleLayer,Identity);
	       break;
      case 3:  IntLin3(c,&E[ind1],&E[ind2],s,t,&Q[g1],R,M,SingleLayer,DoubleLayer);
	       break;
      case 4:  IntLin4(c,&E[ind1],&E[ind2],s,t,&Q[g1],R,M,SingleLayer,DoubleLayer);
               break;
      default: printf("ERROR: CASE != 1,2,3,4 \n");
      }
   }
else 		/* (E[ind1].level > E[ind2].level) */
{  if (dist*(1<<E[ind2].level) >= q)		/* Quadratur */
   {  quadrature_grade_pwl(&g1,&g2,E[ind1].level,E[ind2].level,dist,prec);
      IntLin1(c,&E[ind1],&E[ind2],&RW[g1],&Q[g1],&Q[g2],R,M,SingleLayer,DoubleLayer);
      }
   else		/* unterteile */
   {  element_element_interaction_pwl(&a[  0],P,E,ind1,E[ind2].son[0],RW,Q,R,M,prec,SingleLayer,DoubleLayer,Identity);
      element_element_interaction_pwl(&a[ 48],P,E,ind1,E[ind2].son[1],RW,Q,R,M,prec,SingleLayer,DoubleLayer,Identity);
      element_element_interaction_pwl(&a[ 96],P,E,ind1,E[ind2].son[2],RW,Q,R,M,prec,SingleLayer,DoubleLayer,Identity);
      element_element_interaction_pwl(&a[144],P,E,ind1,E[ind2].son[3],RW,Q,R,M,prec,SingleLayer,DoubleLayer,Identity);

      /* bilde die jeweils vier Integrale */
      for (k=0; k<48; k+=4)
      {  m = 0.25*(a[k+2]+a[48+k+3]+a[96+k]+a[144+k+1]);
         e[0] = 0.5*(a[    k+1]+a[ 48+k  ]);
         e[1] = 0.5*(a[ 48+k+2]+a[ 96+k+1]);
         e[2] = 0.5*(a[ 96+k+3]+a[144+k+2]);
         e[3] = 0.5*(a[144+k  ]+a[    k+3]);

         c[k  ] = 0.5*(a[    k  ]+e[3]+e[0]+m);
         c[k+1] = 0.5*(a[ 48+k+1]+e[0]+e[1]+m);
         c[k+2] = 0.5*(a[ 96+k+2]+e[1]+e[2]+m);
         c[k+3] = 0.5*(a[144+k+3]+e[2]+e[3]+m);
         }
      }
   }

/* Abspeichern des Integrals */
if (ind1 >= ind2) set_integral_pwl(&E[ind1].interaction,ind2,c);
else 
{  permutate(a,c);
   set_integral_pwl(&E[ind2].interaction,ind1,a);
   }
return;
}
