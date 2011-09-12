/**************
 *  Trafos.c  *
 **************/


/*==================================================*
 *  Enthaelt alle Hilfsroutinen fuer die Quadratur  *
 *==================================================*/


#include <math.h>
#include "vector2.h"
#include "trafos.h"


vector2 Kappa(a,b,h)
/* Abbildung vom Einheitsquadrat auf das Element */
vector2		a, b;
double		h;
{
vector2		c;
c.x = a.x + h*b.x;
c.y = a.y + h*b.y;
return(c);
}


vector2 Tau(b_x,b_y,CASE)
/* definiert die Drehungen fuer den Duffy-Trick */
double		b_x, b_y;
unsigned int	CASE;
{
switch (CASE)
{  case 1:  return(vector2_make(1-b_y,  b_x));	/* Drehung um pi/2   */
   case 2:  return(vector2_make(1-b_x,1-b_y));	/* Drehung um pi     */
   case 3:  return(vector2_make(  b_y,1-b_x));	/* Drehung um 3*pi/2 */ 
   default: return(vector2_make(  b_x,  b_y));	/* Identitaet        */
   }
}


unsigned int compare(F1,F2,ind1,ind2)
/* Untersucht zwei Patches (F1[0],...,F1[3]) und (F2[0],...,F2[3])
   auf eine gemeinsame Kante (Funktionsergebnis: 3) oder einen 
   gemeinsamen Punkt (Funktionsergebnis: 4), wobei in diesen Faellen 
   die Drehungen ind1 und ind2 entsprechend der Situation richtig gesetzt
   werden. Im Falle keiner Gemeinsamkeiten ist das Funktionsergebnis 1. */
unsigned int	*F1, *F2;	/* Arrays mit den Eckpunkte der gegebenen Patches */
unsigned int	*ind1, *ind2;	/* Indizes der Drehungen fuer die Patches         */
{
/* zuerst auf einen gemeinsamen Punkt untersuchen */
for (*ind1=0; *ind1<4; (*ind1)++)
{  for (*ind2=0; *ind2<4; (*ind2)++)
   {  if (F1[*ind1] == F2[*ind2])
      {  /* gemeinsamer Punkt -> untersuche auf gemneinsame Kante */
         if (F1[3] == F2[(*ind2+1)%4])			/* dies ist der Sonderfall, denn falls der erste       */
         {  *ind1 = 3;					/* gemeinsame Punkt fuer ind1 = 0 auftritt,            */
            return(3);					/* kann der zweite auch bei ind3 = 3 liegen!           */
            }
         else if (F1[(*ind1+1)%4] == F2[(*ind2+3)%4])  	/* normaler Fall: zweiter Punkt bei ind1+1 bzw. ind2-1 */
         {  *ind2 = (*ind2+3)%4;
            return(3);
            }
         else return(4);
         }
      }
   }
return(1);
}


void quadrature_grade_smooth(g1,g2,m1,m2,dist,prec)
/* bestimmt fuer Elemente der Level m1 und m2 
   mit dem Abstand dist die benoetigten Quadraturgrade
   -> beide Elemente leben auf dem gleichem Patch */
signed int	*g1, *g2;		/* benoetigter Quadraturgrad      */
unsigned int 	m1, m2;			/* Level der zwei Elemente        */
double		dist;			/* Abstand zwischen den Elementen */
double		prec;			/* Quadraturgenauigkeit           */
{
/* Initialisierung */
prec += m1+m2;
dist = (dist < 1) ? log(dist)/log(2) : 0;

/* Quadraturgrad g1 */
*g1 = (signed int) (-0.5*(prec+dist)/(dist+m1+2));
if (*g1 < 0) *g1 = 0;

/* Quadraturgrad g2 */
*g2 = (signed int) (-0.5*(prec+dist)/(dist+m2+2));
if (*g2 < 0) *g2 = 0;

return;
}


void quadrature_grade_nonsmooth(g1,g2,m1,m2,dist,prec)
/* bestimmt fuer Elemente der Level m1 und m2 
   mit dem Abstand dist die benoetigten Quadraturgrade
   -> beide Elemente leben auf verschiedenen Patches */
signed int	*g1, *g2;		/* benoetigter Quadraturgrad      */
unsigned int 	m1, m2;			/* Level der zwei Elemente        */
double		dist;			/* Abstand zwischen den Elementen */
double		prec;			/* Quadraturgenauigkeit           */
{
/* Initialisierung */
prec += m1+m2;
dist = (dist < 1) ? log(dist)/log(2) : 0;

/* Quadraturgrad g1 */
*g1 = (signed int) (-0.5*(prec+2*dist)/(dist+m1+2));
if (*g1 < 0) *g1 = 0;

/* Quadraturgrad g2 */
*g2 = (signed int) (-0.5*(prec+2*dist)/(dist+m2+2));
if (*g2 < 0) *g2 = 0;

return;
}
