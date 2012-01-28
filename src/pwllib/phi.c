/***********
 *  Phi.c  *
 ***********/

 
/*============================================*
 *  Definiert die vier stueckweise linearen   *
 *  Basisfunktionen auf dem Einheitsquadrat.  *
 *  Es ist Phi0([0,0]) = 1,		      *
 *         Phi1([1,0]) = 1, 		      *
 * 	   Phi2([1,1]) = 1,		      *
 *         Phi3([0,1]) = 1.		      *
 *============================================*/


#include <stdlib.h>
#include "vector2.h"
#include "vector3.h"
#include "phi.h"


/* Ansatzfunktion 0 */
double Phi0(a)
vector2	a;
{  return((1-a.x)*(1-a.y));  }


/* Ansatzfunktion 1 */
double Phi1(a)
vector2	a;
{  return(a.x*(1-a.y));  }


/* Ansatzfunktion 2 */
double Phi2(a)
vector2	a;
{  return(a.x*a.y);  }


/* Ansatzfunktion 3 */
double Phi3(a)
vector2	a;
{  return((1-a.x)*a.y);  }


void Phi_times_Phi(c,weight,xi,eta)
/* updated c_{i,j} um weight*phi_i(xi)*phi_j(eta) */ 
double		*c, weight;
vector2		xi, eta;
{
double	a0, a1, a2, a3;
double	b0, b1, b2, b3;

a0 = weight*(1-xi.x)*(1-xi.y);
a1 = weight*   xi.x *(1-xi.y);
a2 = weight*   xi.x *   xi.y ;
a3 = weight*(1-xi.x)*   xi.y ;

b0 = (1-eta.x)*(1-eta.y);
b1 =    eta.x *(1-eta.y);
b2 =    eta.x *   eta.y ;
b3 = (1-eta.x)*   eta.y ;

c[ 0] += a0*b0;
c[ 1] += a0*b1;
c[ 2] += a0*b2;
c[ 3] += a0*b3;
c[ 4] += a1*b0;
c[ 5] += a1*b1;
c[ 6] += a1*b2;
c[ 7] += a1*b3;
c[ 8] += a2*b0;
c[ 9] += a2*b1;
c[10] += a2*b2;
c[11] += a2*b3;
c[12] += a3*b0;
c[13] += a3*b1;
c[14] += a3*b2;
c[15] += a3*b3;
return;
}


void Curl_Phi_times_Curl_Phi(c,weight,xi,eta,dChi_dx_pwl_s,dChi_dy_pwl_s,dChi_dx_t,dChi_dy_t)
/* updated c_{i,j} um weight * < curl[phi_i(xi)],curl[phi_j(eta)] > */ 
double		*c, weight;
vector2		xi, eta;
vector3		dChi_dx_pwl_s, dChi_dy_pwl_s, dChi_dx_t, dChi_dy_t;
{
double		a0, a1, a2, a3;
double		b0, b1, b2, b3;
double		c0, c1, c2, c3;

a0 = weight*vector3_skalp(dChi_dy_pwl_s,dChi_dy_t);
a1 = weight*vector3_skalp(dChi_dy_pwl_s,dChi_dx_pwl_t);
a2 = weight*vector3_skalp(dChi_dx_pwl_s,dChi_dy_pwl_t);
a3 = weight*vector3_skalp(dChi_dx_pwl_s,dChi_dx_t);

b0 = (1-eta.y)*a0 - (1-eta.x)*a1;
b1 = (1-eta.y)*a0 +    eta.x *a1;
b2 =    eta.y *a0 -    eta.x *a1;
b3 =    eta.y *a0 + (1-eta.x)*a1;

c0 = (1-eta.y)*a2 - (1-eta.x)*a3;
c1 = (1-eta.y)*a2 +    eta.x *a3;
c2 =    eta.y *a2 -    eta.x *a3;
c3 =    eta.y *a2 + (1-eta.x)*a3;

c[ 0] += +(1-xi.y)*b0 - (1-xi.x)*c0;
c[ 1] += -(1-xi.y)*b1 + (1-xi.x)*c1;
c[ 2] += -(1-xi.y)*b2 + (1-xi.x)*c2;
c[ 3] += +(1-xi.y)*b3 - (1-xi.x)*c3;

c[ 4] += -(1-xi.y)*b0 -    xi.x *c0;
c[ 5] += +(1-xi.y)*b1 +    xi.x *c1;
c[ 6] += +(1-xi.y)*b2 +    xi.x *c2;
c[ 7] += -(1-xi.y)*b3 -    xi.x *c3;

c[ 8] += -   xi.y *b0 +    xi.x *c0;
c[ 9] += +   xi.y *b1 -    xi.x *c1;
c[10] += +   xi.y *b2 -    xi.x *c2;
c[11] += -   xi.y *b3 +    xi.x *c3;

c[12] += +   xi.y *b0 + (1-xi.x)*c0;
c[13] += -   xi.y *b1 - (1-xi.x)*c1;
c[14] += -   xi.y *b2 - (1-xi.x)*c2;
c[15] += +   xi.y *b3 + (1-xi.x)*c3;
return;
}
