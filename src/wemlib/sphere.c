/**************
 *  Sphere.c  *
 **************/


/*========================================*
 *  Parametrisierung der Einheitssphaere  * 
 *  ueber dem Einheitswuerfel		  *
 *========================================*/


#include <math.h>
#include <stdlib.h>
#include "vector2.h"  
#include "vector3.h"  
#include "gamma.h"


/* Patch 0 */
vector3 Chi0(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.x-0.5;
c.y = s.y-0.5;
c.z = 0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x *= r;
c.y *= r;
c.z *= r;
return(c);
}


vector3 dChi0_dx(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.x-0.5;
c.y = s.y-0.5;
c.z = 0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x = r-c.x*(s.x-0.5)*r*r*r;
c.y =  -c.y*(s.x-0.5)*r*r*r;
c.z =  -c.z*(s.x-0.5)*r*r*r;
return(c);
}


vector3 dChi0_dy(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.x-0.5;
c.y = s.y-0.5;
c.z = 0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x =  -c.x*(s.y-0.5)*r*r*r;
c.y = r-c.y*(s.y-0.5)*r*r*r;
c.z =  -c.z*(s.y-0.5)*r*r*r;
return(c);
}


vector3 n_Chi0(s)
vector2		s;
{  
vector3 	c, dc_dx, dc_dy;
double  	r;

c.x = s.x-0.5;
c.y = s.y-0.5;
c.z = 0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
dc_dx.x = r-c.x*(s.x-0.5)*r*r*r;
dc_dx.y =  -c.y*(s.x-0.5)*r*r*r;
dc_dx.z =  -c.z*(s.x-0.5)*r*r*r;
dc_dy.x =  -c.x*(s.y-0.5)*r*r*r;
dc_dy.y = r-c.y*(s.y-0.5)*r*r*r;
dc_dy.z =  -c.z*(s.y-0.5)*r*r*r;
c.x = dc_dx.y*dc_dy.z - dc_dx.z*dc_dy.y;
c.y = dc_dx.z*dc_dy.x - dc_dx.x*dc_dy.z;
c.z = dc_dx.x*dc_dy.y - dc_dx.y*dc_dy.x;
return(c);
}


/* Patch 1 */
vector3 Chi1(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.y-0.5;
c.y = s.x-0.5;
c.z = -0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x *= r;
c.y *= r;
c.z *= r;
return(c);
}


vector3 dChi1_dx(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.y-0.5;
c.y = s.x-0.5;
c.z = -0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x =  -c.x*(s.x-0.5)*r*r*r;
c.y = r-c.y*(s.x-0.5)*r*r*r;
c.z =  -c.z*(s.x-0.5)*r*r*r;
return(c);
}


vector3 dChi1_dy(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.y-0.5;
c.y = s.x-0.5;
c.z = -0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x = r-c.x*(s.y-0.5)*r*r*r;
c.y =  -c.y*(s.y-0.5)*r*r*r;
c.z =  -c.z*(s.y-0.5)*r*r*r;
return(c);
}


vector3 n_Chi1(s)
vector2		s;
{  
vector3 	c, dc_dx, dc_dy;
double  	r;

c.x = s.y-0.5;
c.y = s.x-0.5;
c.z = -0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
dc_dx.x =  -c.x*(s.x-0.5)*r*r*r;
dc_dx.y = r-c.y*(s.x-0.5)*r*r*r;
dc_dx.z =  -c.z*(s.x-0.5)*r*r*r;
dc_dy.x = r-c.x*(s.y-0.5)*r*r*r;
dc_dy.y =  -c.y*(s.y-0.5)*r*r*r;
dc_dy.z =  -c.z*(s.y-0.5)*r*r*r;
c.x = dc_dx.y*dc_dy.z - dc_dx.z*dc_dy.y;
c.y = dc_dx.z*dc_dy.x - dc_dx.x*dc_dy.z;
c.z = dc_dx.x*dc_dy.y - dc_dx.y*dc_dy.x;
return(c);
}


/* Patch 2 */
vector3 Chi2(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = 0.5;
c.y = s.x-0.5;
c.z = s.y-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x *= r;
c.y *= r;
c.z *= r;
return(c);
}


vector3 dChi2_dx(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = 0.5;
c.y = s.x-0.5;
c.z = s.y-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x =  -c.x*(s.x-0.5)*r*r*r;
c.y = r-c.y*(s.x-0.5)*r*r*r;
c.z =  -c.z*(s.x-0.5)*r*r*r;
return(c);
}


vector3 dChi2_dy(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = 0.5;
c.y = s.x-0.5;
c.z = s.y-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x =  -c.x*(s.y-0.5)*r*r*r;
c.y =  -c.y*(s.y-0.5)*r*r*r;
c.z = r-c.z*(s.y-0.5)*r*r*r;
return(c);
}


vector3 n_Chi2(s)
vector2		s;
{  
vector3 	c, dc_dx, dc_dy;
double  	r;

c.x = 0.5;
c.y = s.x-0.5;
c.z = s.y-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
dc_dx.x =  -c.x*(s.x-0.5)*r*r*r;
dc_dx.y = r-c.y*(s.x-0.5)*r*r*r;
dc_dx.z =  -c.z*(s.x-0.5)*r*r*r;
dc_dy.x =  -c.x*(s.y-0.5)*r*r*r;
dc_dy.y =  -c.y*(s.y-0.5)*r*r*r;
dc_dy.z = r-c.z*(s.y-0.5)*r*r*r;
c.x = dc_dx.y*dc_dy.z - dc_dx.z*dc_dy.y;
c.y = dc_dx.z*dc_dy.x - dc_dx.x*dc_dy.z;
c.z = dc_dx.x*dc_dy.y - dc_dx.y*dc_dy.x;
return(c);
}


/* Patch 3 */
vector3 Chi3(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = -0.5;
c.y = s.y-0.5;
c.z = s.x-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x *= r;
c.y *= r;
c.z *= r;
return(c);
}


vector3 dChi3_dx(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = -0.5;
c.y = s.y-0.5;
c.z = s.x-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x =  -c.x*(s.x-0.5)*r*r*r;
c.y =  -c.y*(s.x-0.5)*r*r*r;
c.z = r-c.z*(s.x-0.5)*r*r*r;
return(c);
}


vector3 dChi3_dy(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = -0.5;
c.y = s.y-0.5;
c.z = s.x-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x =  -c.x*(s.y-0.5)*r*r*r;
c.y = r-c.y*(s.y-0.5)*r*r*r;
c.z =  -c.z*(s.y-0.5)*r*r*r;
return(c);
}


vector3 n_Chi3(s)
vector2		s;
{  
vector3 	c, dc_dx, dc_dy;
double  	r;

c.x = -0.5;
c.y = s.y-0.5;
c.z = s.x-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
dc_dx.x =  -c.x*(s.x-0.5)*r*r*r;
dc_dx.y =  -c.y*(s.x-0.5)*r*r*r;
dc_dx.z = r-c.z*(s.x-0.5)*r*r*r;
dc_dy.x =  -c.x*(s.y-0.5)*r*r*r;
dc_dy.y = r-c.y*(s.y-0.5)*r*r*r;
dc_dy.z =  -c.z*(s.y-0.5)*r*r*r;
c.x = dc_dx.y*dc_dy.z - dc_dx.z*dc_dy.y;
c.y = dc_dx.z*dc_dy.x - dc_dx.x*dc_dy.z;
c.z = dc_dx.x*dc_dy.y - dc_dx.y*dc_dy.x;
return(c);
}


/* Patch 4 */
vector3 Chi4(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.y-0.5;
c.y = 0.5;
c.z = s.x-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x *= r;
c.y *= r;
c.z *= r;
return(c);
}


vector3 dChi4_dx(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.y-0.5;
c.y = 0.5;
c.z = s.x-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x =  -c.x*(s.x-0.5)*r*r*r;
c.y =  -c.y*(s.x-0.5)*r*r*r;
c.z = r-c.z*(s.x-0.5)*r*r*r;
return(c);
}


vector3 dChi4_dy(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.y-0.5;
c.y = 0.5;
c.z = s.x-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x = r-c.x*(s.y-0.5)*r*r*r;
c.y =  -c.y*(s.y-0.5)*r*r*r;
c.z =  -c.z*(s.y-0.5)*r*r*r;
return(c);
}


vector3 n_Chi4(s)
vector2		s;
{  
vector3 	c, dc_dx, dc_dy;
double  	r;

c.x = s.y-0.5;
c.y = 0.5;
c.z = s.x-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
dc_dx.x =  -c.x*(s.x-0.5)*r*r*r;
dc_dx.y =  -c.y*(s.x-0.5)*r*r*r;
dc_dx.z = r-c.z*(s.x-0.5)*r*r*r;
dc_dy.x = r-c.x*(s.y-0.5)*r*r*r;
dc_dy.y =  -c.y*(s.y-0.5)*r*r*r;
dc_dy.z =  -c.z*(s.y-0.5)*r*r*r;
c.x = dc_dx.y*dc_dy.z - dc_dx.z*dc_dy.y;
c.y = dc_dx.z*dc_dy.x - dc_dx.x*dc_dy.z;
c.z = dc_dx.x*dc_dy.y - dc_dx.y*dc_dy.x;
return(c);
}


/* Patch 5 */
vector3 Chi5(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.x-0.5;
c.y = -0.5;
c.z = s.y-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x *= r;
c.y *= r;
c.z *= r;
return(c);
}


vector3 dChi5_dx(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.x-0.5;
c.y = -0.5;
c.z = s.y-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x = r-c.x*(s.x-0.5)*r*r*r;
c.y =  -c.y*(s.x-0.5)*r*r*r;
c.z =  -c.z*(s.x-0.5)*r*r*r;
return(c);
}


vector3 dChi5_dy(s)
vector2		s;
{  
vector3 	c;
double  	r;

c.x = s.x-0.5;
c.y = -0.5;
c.z = s.y-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
c.x =  -c.x*(s.y-0.5)*r*r*r;
c.y =  -c.y*(s.y-0.5)*r*r*r;
c.z = r-c.z*(s.y-0.5)*r*r*r;
return(c);
}


vector3 n_Chi5(s)
vector2		s;
{  
vector3 	c, dc_dx, dc_dy;
double  	r;

c.x = s.x-0.5;
c.y = -0.5;
c.z = s.y-0.5;
r = 1/sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
dc_dx.x = r-c.x*(s.x-0.5)*r*r*r;
dc_dx.y =  -c.y*(s.x-0.5)*r*r*r;
dc_dx.z =  -c.z*(s.x-0.5)*r*r*r;
dc_dy.x =  -c.x*(s.y-0.5)*r*r*r;
dc_dy.y =  -c.y*(s.y-0.5)*r*r*r;
dc_dy.z = r-c.z*(s.y-0.5)*r*r*r;
c.x = dc_dx.y*dc_dy.z - dc_dx.z*dc_dy.y;
c.y = dc_dx.z*dc_dy.x - dc_dx.x*dc_dy.z;
c.z = dc_dx.x*dc_dy.y - dc_dx.y*dc_dy.x;
return(c);
}


unsigned int init_p()
/* Initialisierung: Liefert als Funktionsergebnis die 
   Anzahl der Parametergebiete */
{  return(6);  }   
   
   
void init_Chi(Chi)
/* allokiert den noetigen Speicherplatz fuer die 
   Parametrisierung und definiert Chi[0],...,Chi[5] */
parametrix	**Chi;
{  
(*Chi) = (parametrix*) malloc(6*sizeof(parametrix));

(*Chi)[0].f = Chi0;
(*Chi)[0].df_dx = dChi0_dx;
(*Chi)[0].df_dy = dChi0_dy;
(*Chi)[0].n_f = n_Chi0;

(*Chi)[1].f = Chi1;
(*Chi)[1].df_dx = dChi1_dx;
(*Chi)[1].df_dy = dChi1_dy;
(*Chi)[1].n_f = n_Chi1;

(*Chi)[2].f = Chi2;
(*Chi)[2].df_dx = dChi2_dx;
(*Chi)[2].df_dy = dChi2_dy;
(*Chi)[2].n_f = n_Chi2;

(*Chi)[3].f = Chi3;
(*Chi)[3].df_dx = dChi3_dx;
(*Chi)[3].df_dy = dChi3_dy;
(*Chi)[3].n_f = n_Chi3;

(*Chi)[4].f = Chi4;
(*Chi)[4].df_dx = dChi4_dx;
(*Chi)[4].df_dy = dChi4_dy;
(*Chi)[4].n_f = n_Chi4;

(*Chi)[5].f = Chi5;
(*Chi)[5].df_dx = dChi5_dx;
(*Chi)[5].df_dy = dChi5_dy;
(*Chi)[5].n_f = n_Chi5;

return;
}


void free_Chi(Chi)
/* gibt den Speicherplatz fuer die Parametrisierung frei */
parametrix	**Chi;
{  free(*Chi);  }
