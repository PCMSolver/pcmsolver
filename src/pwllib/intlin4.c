/***************
 *  IntLin4.c  *
 ***************/


#include <math.h>
#include <string.h>
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "trafos.h"
#include "kern.h"
#include "phi.h"
#include "interpolate_pwl.h"
#include "cubature.h"
#include "intlin4.h"


void IntLin4(c,element1,element2,ind_s,ind_t,Q,P,M,SingleLayer,DoubleLayer)
/* GEMEINSAME ECKE IM NULLPUNKT -> modifiziertes Skalarprodukt */
double		*c;
element		*element1, *element2;
unsigned int	ind_s, ind_t;
cubature	*Q;
vector3		****P;
unsigned int	M;
double		SingleLayer(), DoubleLayer();
{
unsigned int	i, j;
double		d1, d2, d3, w;
vector2		s, t, xi, eta, a, u, a1, a2, b1, b2;
vector3		x1, n_x1, x2, n_x2, y1, n_y1, y2, n_y2, z, n_z;
double		h = 1./(1 << element1->level);

memset(c,0,48*sizeof(double));
s = vector2_make(h*element1->index_s,h*element1->index_t);
t = vector2_make(h*element2->index_s,h*element2->index_t);

for (i=0; i<Q->nop; i++)
{  xi = Q->xi[i];
   w  = h * h * pow(xi.x,3) * Q->w[i];
   xi.y *= xi.x;
   a1 = Tau(xi.x,xi.y,ind_s);
   a2 = Tau(xi.y,xi.x,ind_s);
   b1 = Tau(xi.x,xi.y,ind_t);
   b2 = Tau(xi.y,xi.x,ind_t);
   
   u = Kappa(s,a1,h);
   x1 = Chi_pwl(u,P[element1->patch],M);
   n_x1 = n_Chi_pwl(u,P[element1->patch],M);
   
   u = Kappa(s,a2,h);
   x2 = Chi_pwl(u,P[element1->patch],M);
   n_x2 = n_Chi_pwl(u,P[element1->patch],M);
      
   u = Kappa(t,b1,h);
   y1 = Chi_pwl(u,P[element2->patch],M);
   n_y1 = n_Chi_pwl(u,P[element2->patch],M);
   
   u = Kappa(t,b2,h);
   y2 = Chi_pwl(u,P[element2->patch],M);
   n_y2 = n_Chi_pwl(u,P[element2->patch],M);

   for (j=0; j<Q->nop; j++)  
   {  eta = vector2_Smul(xi.x,Q->xi[j]);
            
      a = Tau(eta.x,eta.y,ind_t);
      u = Kappa(t,a,h);
      z = Chi_pwl(u,P[element2->patch],M);
      n_z = n_Chi_pwl(u,P[element2->patch],M);

      d1 = w * Q->w[j] * SingleLayer(x1,z);
      d2 = w * Q->w[j] * DoubleLayer(x1,z,n_z );
      d3 = w * Q->w[j] * DoubleLayer(z,x1,n_x1);
      Phi_times_Phi(&c[ 0],d1,a1,a);
      Phi_times_Phi(&c[16],d2,a1,a);
      Phi_times_Phi(&c[32],d3,a1,a);
      
      d1 = w * Q->w[j] * SingleLayer(x2,z);
      d2 = w * Q->w[j] * DoubleLayer(x2,z,n_z );
      d3 = w * Q->w[j] * DoubleLayer(z,x2,n_x2);
      Phi_times_Phi(&c[ 0],d1,a2,a);
      Phi_times_Phi(&c[16],d2,a2,a);
      Phi_times_Phi(&c[32],d3,a2,a);

      a = Tau(eta.x,eta.y,ind_s);
      u = Kappa(s,a,h);
      z = Chi_pwl(u,P[element1->patch],M);
      n_z = n_Chi_pwl(u,P[element1->patch],M);
      
      d1 = w * Q->w[j] * SingleLayer(z,y1);
      d2 = w * Q->w[j] * DoubleLayer(z,y1,n_y1);
      d3 = w * Q->w[j] * DoubleLayer(y1,z,n_z );
      Phi_times_Phi(&c[ 0],d1,a,b1);
      Phi_times_Phi(&c[16],d2,a,b1);
      Phi_times_Phi(&c[32],d3,a,b1);

      d1 = w * Q->w[j] * SingleLayer(z,y2);
      d2 = w * Q->w[j] * DoubleLayer(z,y2,n_y2);
      d3 = w * Q->w[j] * DoubleLayer(y2,z,n_z );
      Phi_times_Phi(&c[ 0],d1,a,b2);
      Phi_times_Phi(&c[16],d2,a,b2);
      Phi_times_Phi(&c[32],d3,a,b2);
      } 
   } 	 
return;
}
