/**************
 *  Trafos.c  *
 **************/


/*==================================================*
 *  Enthaelt alle Hilfsroutinen fuer die Quadratur  *
 *==================================================*/


#include <math.h>
#include "vector2.h"
#include "vector3.h"
#include "constants.h"
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


unsigned int compare(P,F1,F2,ind1,ind2)
/* Untersucht zwei Patches (F1[0],...,F1[3]) und (F2[0],...,F2[3])
   auf eine gemeinsame Kante (Funktionsergebnis: 3) oder einen 
   gemeinsamen Punkt (Funktionsergebnis: 4), wobei in diesen Faellen 
   die Drehungen ind1 und ind2 entsprechend der Situation richtig gesetzt
   werden. Im Falle keiner Gemeinsamkeiten ist das Funktionsergebnis 1. */

vector3		*P;		/* Punkteliste					  */
unsigned int	*F1, *F2;	/* Arrays mit den Eckpunkte der gegebenen Patches */
unsigned int	*ind1, *ind2;	/* Indizes der Drehungen fuer die Patches         */
{
vector3		d;		/* Differenzvektor der Elementeckpunkte           */

/* zuerst auf einen gemeinsamen Punkt untersuchen */
for (*ind1=0; *ind1<4; (*ind1)++)
{  for (*ind2=0; *ind2<4; (*ind2)++)
   {  d.x = P[F1[*ind1]].x-P[F2[*ind2]].x;
      d.y = P[F1[*ind1]].y-P[F2[*ind2]].y;
      d.z = P[F1[*ind1]].z-P[F2[*ind2]].z;
      if (d.x*d.x+d.y*d.y+d.z*d.z < tol)
      {  /* gemeinsamer Punkt -> untersuche auf gemneinsame Kante */
	 d.x = P[F1[(*ind1+1)%4]].x-P[F2[(*ind2+3)%4]].x;
         d.y = P[F1[(*ind1+1)%4]].y-P[F2[(*ind2+3)%4]].y;
         d.z = P[F1[(*ind1+1)%4]].z-P[F2[(*ind2+3)%4]].z;
	 if (d.x*d.x+d.y*d.y+d.z*d.z < tol)
         {  *ind2 = (*ind2+3)%4;	/* normaler Fall: zweiter Punkt bei ind1+1 bzw. ind2-1 */
            return(3);
            }
         else if (*ind1 == 0)
	 {  d.x = P[F1[3]].x-P[F2[(*ind2+1)%4]].x;	/* dies ist der Sonderfall, denn falls der erste gemeinsame Punkt    */ 
            d.y = P[F1[3]].y-P[F2[(*ind2+1)%4]].y;	/* fuer ind1 = 0 auftritt, kann der zweite auch bei ind3 = 3 liegen! */	
            d.z = P[F1[3]].z-P[F2[(*ind2+1)%4]].z;
            if (d.x*d.x+d.y*d.y+d.z*d.z < tol)
            {  *ind1 = 3;
               return(3);
               }
            }
         return(4);
         }
      }
   }
return(1);
}


void quadrature_grade(g1,g2,m1,m2,dist,prec)
/* bestimmt fuer Elemente der Level m1 und m2 
   mit dem Abstand dist die benoetigten Quadraturgrade */
signed int	*g1, *g2;		/* benoetigter Quadraturgrad      */
unsigned int 	m1, m2;			/* Level der zwei Elemente        */
double		dist;			/* Abstand zwischen den Elementen */
double		prec;			/* Quadraturgenauigkeit           */
{
/* Initialisierung */
dist = (dist < 1) ? log(dist)/log(2) : 0;

/* Quadraturgrad g1 */
*g1 = (signed int) (-0.5*(prec+m2)/(dist+m1+2));
if (*g1 < 1) *g1 = 1;

/* Quadraturgrad g2 */
*g2 = (signed int) (-0.5*(prec+m1)/(dist+m2+2));
if (*g2 < 1) *g2 = 1;

return;
}
