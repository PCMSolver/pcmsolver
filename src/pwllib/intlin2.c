/***************
 *  IntLin2.c  *
 ***************/


#include <string.h>
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "kern.h"
#include "phi.h"
#include "cubature.h"
#include "interpolate_pwl.h"
#include "intlin2.h"


void IntLin2(c,element1,Q,P,M,SingleLayer,DoubleLayer,Identity)
/* GLEICHE PATCHES [0,1]^2 -> modifiziertes Skalarprodukt */
double		*c;
element		*element1;
cubature	*Q;
vector3		****P;
unsigned int	M;
double		SingleLayer(), DoubleLayer();
double		Identity;
{
unsigned int	i, j;
double		d1, d2, d3, w;
double		t1, t2, t3, t4;
vector2		s, xi, eta, a, b;
vector3		x, y;
double		h = 1./(1 << element1->level);

memset(c,0,32*sizeof(double));
s = vector2_make(h*element1->index_s,h*element1->index_t);

for (i=0; i<Q->nop; i++)
{  xi = Q->xi[i];
   w  = h * h * Q->w[i] * xi.x * (1-xi.x) * (1-xi.x*xi.y);
   for (j=0; j<Q->nop; j++)
   {  eta = Q->xi[j];
      t1 = eta.x*(1-xi.x);
      t2 = eta.y*(1-xi.x*xi.y);
      t3 = t1+xi.x;
      t4 = t2+xi.x*xi.y;

      a.x = s.x + h*t1;
      a.y = s.y + h*t2;
      b.x = s.x + h*t3;
      b.y = s.y + h*t4;
      x = Chi_pwl(a,P[element1->patch],M);
      y = Chi_pwl(b,P[element1->patch],M);
      d1 = w * Q->w[j] * SingleLayer(x,y);
      d2 = w * Q->w[j] * DoubleLayer(x,y,n_Chi_pwl(b,P[element1->patch],M));
      d3 = w * Q->w[j] * DoubleLayer(y,x,n_Chi_pwl(a,P[element1->patch],M));
      Phi_times_Phi(&c[ 0],d1,vector2_make(t1,t2),vector2_make(t3,t4));
      Phi_times_Phi(&c[ 0],d1,vector2_make(t3,t4),vector2_make(t1,t2));
      Phi_times_Phi(&c[16],d2,vector2_make(t1,t2),vector2_make(t3,t4));
      Phi_times_Phi(&c[16],d3,vector2_make(t3,t4),vector2_make(t1,t2));

      a.y = s.y + h*t4;
      b.y = s.y + h*t2;
      x = Chi_pwl(a,P[element1->patch],M);
      y = Chi_pwl(b,P[element1->patch],M);
      d1 = w * Q->w[j] * SingleLayer(x,y);
      d2 = w * Q->w[j] * DoubleLayer(x,y,n_Chi_pwl(b,P[element1->patch],M));
      d3 = w * Q->w[j] * DoubleLayer(y,x,n_Chi_pwl(a,P[element1->patch],M));
      Phi_times_Phi(&c[ 0],d1,vector2_make(t1,t4),vector2_make(t3,t2));
      Phi_times_Phi(&c[ 0],d1,vector2_make(t3,t2),vector2_make(t1,t4));
      Phi_times_Phi(&c[16],d2,vector2_make(t1,t4),vector2_make(t3,t2));
      Phi_times_Phi(&c[16],d3,vector2_make(t3,t2),vector2_make(t1,t4));

      a.x = s.x + h*t2;
      a.y = s.y + h*t1;
      b.x = s.x + h*t4;
      b.y = s.y + h*t3;
      x = Chi_pwl(a,P[element1->patch],M);
      y = Chi_pwl(b,P[element1->patch],M);
      d1 = w * Q->w[j] * SingleLayer(x,y);
      d2 = w * Q->w[j] * DoubleLayer(x,y,n_Chi_pwl(b,P[element1->patch],M));
      d3 = w * Q->w[j] * DoubleLayer(y,x,n_Chi_pwl(a,P[element1->patch],M));
      Phi_times_Phi(&c[ 0],d1,vector2_make(t2,t1),vector2_make(t4,t3));
      Phi_times_Phi(&c[ 0],d1,vector2_make(t4,t3),vector2_make(t2,t1));
      Phi_times_Phi(&c[16],d2,vector2_make(t2,t1),vector2_make(t4,t3));
      Phi_times_Phi(&c[16],d3,vector2_make(t4,t3),vector2_make(t2,t1));

      a.y = s.y + h*t3;
      b.y = s.y + h*t1;
      x = Chi_pwl(a,P[element1->patch],M);
      y = Chi_pwl(b,P[element1->patch],M);
      d1 = w * Q->w[j] * SingleLayer(x,y);
      d2 = w * Q->w[j] * DoubleLayer(x,y,n_Chi_pwl(b,P[element1->patch],M));
      d3 = w * Q->w[j] * DoubleLayer(y,x,n_Chi_pwl(a,P[element1->patch],M));
      Phi_times_Phi(&c[ 0],d1,vector2_make(t2,t3),vector2_make(t4,t1));
      Phi_times_Phi(&c[ 0],d1,vector2_make(t4,t1),vector2_make(t2,t3));
      Phi_times_Phi(&c[16],d2,vector2_make(t2,t3),vector2_make(t4,t1));
      Phi_times_Phi(&c[16],d3,vector2_make(t4,t1),vector2_make(t2,t3));
      }
   }

/* bilde + Identity */
w = Identity;
c[16] += w/9;
c[17] += w/18;
c[18] += w/36;
c[19] += w/18;
c[20] += w/18;
c[21] += w/9;
c[22] += w/18; 
c[23] += w/36;
c[24] += w/36;
c[25] += w/18;
c[26] += w/9;
c[27] += w/18;
c[28] += w/18;
c[29] += w/36;
c[30] += w/18;
c[31] += w/9;

/* transponierter Eintrag */
c[32] = c[16];
c[33] = c[20];
c[34] = c[24];
c[35] = c[28];
c[36] = c[17];
c[37] = c[21];
c[38] = c[25];
c[39] = c[29];
c[40] = c[18];
c[41] = c[22];
c[42] = c[26];
c[43] = c[30];
c[44] = c[19];
c[45] = c[23];
c[46] = c[27];
c[47] = c[31];
return;
}
