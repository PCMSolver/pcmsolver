/***************
 *  IntKon3.c  *
 ***************/


#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "trafos.h"
#include "kern.h"
#include "cubature.h"
#include "interpolate.h"
#include "intkon3.h"


void IntKon3(c,element1,element2,ind_s,ind_t,Q,P,M,SingleLayer,DoubleLayer)
/* GEMEINSAME KANTE [0,1] -> modifiziertes Skalarprodukt */
double		*c;
element		*element1, *element2;
unsigned int	ind_s, ind_t;
cubature	*Q;
vector3		****P;
unsigned int	M;
double		SingleLayer(), DoubleLayer();
{
unsigned int	i, j;
double		w, d11, d12, d21, d22, d31, d32, t1, t2, t3, t4;
vector2		s, t, xi, eta, a, b;
vector3		x, y;
double		h = 1./(1 << element1->level);

c[0] = c[1] = c[2] = 0;
s = vector2_make(h*element1->index_s,h*element1->index_t);
t = vector2_make(h*element2->index_s,h*element2->index_t);
for (i=0; i<Q->nop; i++)
{  xi = Q->xi[i];
   w  = xi.y * xi.y * Q->w[i];
   t1 = xi.x*(1-xi.y);
   t2 = (1-xi.x)*(1-xi.y);
   
   for (j=0; j<Q->nop; j++)  
   {  eta = vector2_Smul(xi.y,Q->xi[j]);
      t3 = xi.x*(1-eta.x);
      t4 = (1-xi.x)*(1-eta.x);
      
      a = Kappa(s,Tau(t1,eta.x,ind_s),h);
      b = Kappa(t,Tau(t2,eta.y,ind_t),h);
      x = Chi(a,P[element1->patch],M);
      y = Chi(b,P[element2->patch],M);
      d11 = SingleLayer(x,y);
      d21 = DoubleLayer(x,y,n_Chi(b,P[element2->patch],M));
      d31 = DoubleLayer(y,x,n_Chi(a,P[element1->patch],M));

      a = Kappa(s,Tau(1-t1,eta.x,ind_s),h);
      b = Kappa(t,Tau(1-t2,eta.y,ind_t),h);
      x = Chi(a,P[element1->patch],M);
      y = Chi(b,P[element2->patch],M);
      d11 += SingleLayer(x,y);
      d21 += DoubleLayer(x,y,n_Chi(b,P[element2->patch],M));
      d31 += DoubleLayer(y,x,n_Chi(a,P[element1->patch],M));

      a = Kappa(s,Tau(t3, xi.y,ind_s),h);
      b = Kappa(t,Tau(t4,eta.y,ind_t),h);
      x = Chi(a,P[element1->patch],M);
      y = Chi(b,P[element2->patch],M);
      d12 = SingleLayer(x,y);
      d22 = DoubleLayer(x,y,n_Chi(b,P[element2->patch],M));
      d32 = DoubleLayer(y,x,n_Chi(a,P[element1->patch],M));

      a = Kappa(s,Tau(1-t3, xi.y,ind_s),h);
      b = Kappa(t,Tau(1-t4,eta.y,ind_t),h);
      x = Chi(a,P[element1->patch],M);
      y = Chi(b,P[element2->patch],M);
      d12 += SingleLayer(x,y);
      d22 += DoubleLayer(x,y,n_Chi(b,P[element2->patch],M));
      d32 += DoubleLayer(y,x,n_Chi(a,P[element1->patch],M));

      a = Kappa(s,Tau(t4,eta.y,ind_s),h);
      b = Kappa(t,Tau(t3, xi.y,ind_t),h);
      x = Chi(a,P[element1->patch],M);
      y = Chi(b,P[element2->patch],M);
      d12 += SingleLayer(x,y);
      d22 += DoubleLayer(x,y,n_Chi(b,P[element2->patch],M));
      d32 += DoubleLayer(y,x,n_Chi(a,P[element1->patch],M));

      a = Kappa(s,Tau(1-t4,eta.y,ind_s),h);
      b = Kappa(t,Tau(1-t3, xi.y,ind_t),h);
      x = Chi(a,P[element1->patch],M);
      y = Chi(b,P[element2->patch],M);
      d12 += SingleLayer(x,y);
      d22 += DoubleLayer(x,y,n_Chi(b,P[element2->patch],M));
      d32 += DoubleLayer(y,x,n_Chi(a,P[element1->patch],M));
      
      c[0] += w * Q->w[j] * ((1-xi.y)*d11 + (1-eta.x)*d12);
      c[1] += w * Q->w[j] * ((1-xi.y)*d21 + (1-eta.x)*d22);
      c[2] += w * Q->w[j] * ((1-xi.y)*d31 + (1-eta.x)*d32);
      } 
   } 	 

c[0] *= h*h;		/* L^2-normiert */
c[1] *= h*h;
c[2] *= h*h;
return;
}
