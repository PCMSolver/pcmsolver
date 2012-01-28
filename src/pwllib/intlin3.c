/***************
 *  IntLin3.c  *
 ***************/


#include <string.h>
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "trafos.h"
#include "kern.h"
#include "phi.h"
#include "interpolate.h"
#include "cubature.h"
#include "intlin3.h"


void IntLin3(c,element1,element2,ind_s,ind_t,Q,P,M,SingleLayer,DoubleLayer)
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
double		d1, d2, d3, w, t1, t2, t3, t4;
vector2		s, t, xi, eta, a, b, u, v;
vector3		x, y;
double		h = 1./(1 << element1->level);

memset(c,0,48*sizeof(double));
s = vector2_make(h*element1->index_s,h*element1->index_t);
t = vector2_make(h*element2->index_s,h*element2->index_t);

for (i=0; i<Q->nop; i++)
{  xi = Q->xi[i];
   w  = h * h * xi.y * xi.y * Q->w[i];
   t1 = xi.x*(1-xi.y);
   t2 = (1-xi.x)*(1-xi.y);
   
   for (j=0; j<Q->nop; j++)  
   {  eta = vector2_Smul(xi.y,Q->xi[j]);
      t3 = xi.x*(1-eta.x);
      t4 = (1-xi.x)*(1-eta.x);

      a = Tau(t1,eta.x,ind_s);
      b = Tau(t2,eta.y,ind_t);
      u = Kappa(s,a,h);
      v = Kappa(t,b,h);
      x = Chi(u,P[element1->patch],M);
      y = Chi(v,P[element2->patch],M);
      d1 = w * Q->w[j] * (1-xi.y) * SingleLayer(x,y);
      d2 = w * Q->w[j] * (1-xi.y) * DoubleLayer(x,y,n_Chi(v,P[element2->patch],M));
      d3 = w * Q->w[j] * (1-xi.y) * DoubleLayer(y,x,n_Chi(u,P[element1->patch],M));
      Phi_times_Phi(&c[ 0],d1,a,b);
      Phi_times_Phi(&c[16],d2,a,b);
      Phi_times_Phi(&c[32],d3,a,b);

      a = Tau(1-t1,eta.x,ind_s);
      b = Tau(1-t2,eta.y,ind_t);
      u = Kappa(s,a,h);
      v = Kappa(t,b,h);
      x = Chi(u,P[element1->patch],M);
      y = Chi(v,P[element2->patch],M);
      d1 = w * Q->w[j] * (1-xi.y) * SingleLayer(x,y);
      d2 = w * Q->w[j] * (1-xi.y) * DoubleLayer(x,y,n_Chi(v,P[element2->patch],M));
      d3 = w * Q->w[j] * (1-xi.y) * DoubleLayer(y,x,n_Chi(u,P[element1->patch],M));
      Phi_times_Phi(&c[ 0],d1,a,b);
      Phi_times_Phi(&c[16],d2,a,b);
      Phi_times_Phi(&c[32],d3,a,b);

      a = Tau(t3, xi.y,ind_s);
      b = Tau(t4,eta.y,ind_t);
      u = Kappa(s,a,h);
      v = Kappa(t,b,h);
      x = Chi(u,P[element1->patch],M);
      y = Chi(v,P[element2->patch],M);
      d1 = w * Q->w[j] * (1-eta.x) * SingleLayer(x,y);
      d2 = w * Q->w[j] * (1-eta.x) * DoubleLayer(x,y,n_Chi(v,P[element2->patch],M));
      d3 = w * Q->w[j] * (1-eta.x) * DoubleLayer(y,x,n_Chi(u,P[element1->patch],M));
      Phi_times_Phi(&c[ 0],d1,a,b);
      Phi_times_Phi(&c[16],d2,a,b);
      Phi_times_Phi(&c[32],d3,a,b);

      a = Tau(1-t3, xi.y,ind_s);
      b = Tau(1-t4,eta.y,ind_t);
      u = Kappa(s,a,h);
      v = Kappa(t,b,h);
      x = Chi(u,P[element1->patch],M);
      y = Chi(v,P[element2->patch],M);
      d1 = w * Q->w[j] * (1-eta.x) * SingleLayer(x,y);
      d2 = w * Q->w[j] * (1-eta.x) * DoubleLayer(x,y,n_Chi(v,P[element2->patch],M));
      d3 = w * Q->w[j] * (1-eta.x) * DoubleLayer(y,x,n_Chi(u,P[element1->patch],M));
      Phi_times_Phi(&c[ 0],d1,a,b);
      Phi_times_Phi(&c[16],d2,a,b);
      Phi_times_Phi(&c[32],d3,a,b);

      a = Tau(t4,eta.y,ind_s);
      b = Tau(t3, xi.y,ind_t);
      u = Kappa(s,a,h);
      v = Kappa(t,b,h);
      x = Chi(u,P[element1->patch],M);
      y = Chi(v,P[element2->patch],M);
      d1 = w * Q->w[j] * (1-eta.x) * SingleLayer(x,y);
      d2 = w * Q->w[j] * (1-eta.x) * DoubleLayer(x,y,n_Chi(v,P[element2->patch],M));
      d3 = w * Q->w[j] * (1-eta.x) * DoubleLayer(y,x,n_Chi(u,P[element1->patch],M));
      Phi_times_Phi(&c[ 0],d1,a,b);
      Phi_times_Phi(&c[16],d2,a,b);
      Phi_times_Phi(&c[32],d3,a,b);

      a = Tau(1-t4,eta.y,ind_s);
      b = Tau(1-t3, xi.y,ind_t);
      u = Kappa(s,a,h);
      v = Kappa(t,b,h);
      x = Chi(u,P[element1->patch],M);
      y = Chi(v,P[element2->patch],M);
      d1 = w * Q->w[j] * (1-eta.x) * SingleLayer(x,y);
      d2 = w * Q->w[j] * (1-eta.x) * DoubleLayer(x,y,n_Chi(v,P[element2->patch],M));
      d3 = w * Q->w[j] * (1-eta.x) * DoubleLayer(y,x,n_Chi(u,P[element1->patch],M));
      Phi_times_Phi(&c[ 0],d1,a,b);
      Phi_times_Phi(&c[16],d2,a,b);
      Phi_times_Phi(&c[32],d3,a,b);
      } 
   } 	 
return;
}
