/************
 *  data.c  *
 ************/

 
#include <math.h>
#include "vector3.h"
#include "data.h"


/* number of charges */
const unsigned int k = 12;

/* position of charges: vectors (x,y,z) */
const vector3 x[12] = {
  { 5.274,   1.999,  -8.568 },
  { 6.627,   2.018,  -8.209 },
  { 7.366,   0.829,  -8.202 },
  { 6.752,  -0.379,  -8.554 },
  { 5.399,  -0.398,  -8.912 },
  { 4.660,   0.791,  -8.919 },
  { 4.704,   2.916,  -8.573 },
  { 7.101,   2.950,  -7.938 },
  { 8.410,   0.844,  -7.926 },
  { 7.322,  -1.296,  -8.548 },
  { 4.925,  -1.330,  -9.183 },
  { 3.616,   0.776,  -9.196 }
};

/* charges: double values */
const double alpha[12]={6,6,6, 6,6,6, 1,1,1, 1,1,1}; 


double f(a)
vector3		a;
{
double		c = 0;
unsigned int	i;
for (i=0; i<k; i++) c += alpha[i]/vector3_norm(vector3_make(a.x-x[i].x,a.y-x[i].y,a.z-x[i].z));
return(c);
}


vector3 df(a)
vector3         a;
{
unsigned int	i;
vector3		c, r, v;
c.x = c.y = c.z = 0;
for (i=0; i<k; i++)
{  r = vector3_make(a.x-x[i].x,a.y-x[i].y,a.z-x[i].z);
   v = vector3_Smul(alpha[i]/pow(vector3_norm(r),3),r);
   c.x -= v.x;
   c.y -= v.y;
   c.z -= v.z;
   }
return(c);
}
