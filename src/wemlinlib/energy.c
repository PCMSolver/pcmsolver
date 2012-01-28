/**************
 *  Energy.c  *
 **************/


/*=====================================*
 * Berechnet die potentielle Energie.  *
 *=====================================*/


#include <math.h>
#include <stdio.h>
#include "vector2.h"
#include "vector3.h"
#include "data.h"
#include "cubature.h"
#include "interpolate.h"
#include "gauss_square.h"
#include "energy.h"
#include "kern.h"
#include "phi.h"


#if !defined pi
  #define pi 3.1415926535897932385
#endif  


double energy(u,F,T,p,m)
double		*u; 		/* vorgegebene Dichte			       */
unsigned int	**F;		/* Patchliste				       */
vector3		****T;		/* Koeffizienten zur Oberflaecheninterpolation */
unsigned int	p;	 	/* Anzahl der Parametergebiete                 */
unsigned int	m;		/* Zahl der Level                              */
{
unsigned int	n = 1 << m;	/* n*n Patches pro Parametergebiet             */
unsigned int	i1, i2, i3;	/* Laufindizes fuer Ansatzfunktion   	       */
unsigned int	zi; 		/* Zeilenindex hieraus: zi = i1*(n*n)+i2*n+i3  */
cubature        *Q;		/* Kubatur-Formeln                             */
unsigned int	g = 0;		/* Quadraturgrad                               */
unsigned int	l;	 	/* Laufindex fuer Quadratur		       */
vector2		s;		/* Linker, unterer Eckpunkt des Patches zi     */
vector2		t;		/* Auswertepunkte der Gaussquadratur	       */
double		U;		/* Auswertung der Dichte im Quadraturpunkt     */
double		E;		/* Energie                                     */
double 		h = 1./n;	/* Schrittweite   			       */
double		C = 0;		/* surface meassure                            */
double		D = 0;		/* total charge                                */

/* Initialisierung */
init_Gauss_Square(&Q,g+1);	/* Kubatur-Formeln */

/* Berechnung des Fehlers */
E = zi = 0;
for (i1=0; i1<p; i1++)
{  s.y = 0;
   for (i2=0; i2<n; i2++)
   {  s.x = 0;
      for (i3=0; i3<n; i3++) 	/* zeilenweise Durchnumerierung der Patches zi = (i1,i2,i3) */
      {  for (l=0; l<Q[g].nop; l++)
         {  t = vector2_add(s,vector2_Smul(h,Q[g].xi[l]));
	    U = u[F[zi][0]] * Phi0(Q[g].xi[l]) \
              + u[F[zi][1]] * Phi1(Q[g].xi[l]) \
	      + u[F[zi][2]] * Phi2(Q[g].xi[l]) \
	      + u[F[zi][3]] * Phi3(Q[g].xi[l]) ;
	    E += Q[g].w[l]*U*f(Chi(t,T[i1],m));
	    C += Q[g].w[l]*U;
	    }
	 s.x += h;
	 zi++;
	 }
      s.y += h;
      }
   }
E = -0.5*h*E;	/* correct scaling */

/* Datenausgabe */
printf("Computing the energy:            %g\n",E);

for (l=0; l<k; l++) D += alpha[l];	/* total charge */
/* C = total charge * (epsilon-1)/epsilon */
printf("Testing total charge (= 0):      %g\n",h*C-D*(epsilon-1)/epsilon);
free_Gauss_Square(&Q,g+1);
return(E);
}
