/*******************
 *  Interpolate.c  *
 *******************/


/*=============================================================*
 *  Enthaelt alle Routinen zur Interpolation der Oberflaeche.  *
 *=============================================================*/


#include <stdlib.h>
#include "vector2.h"
#include "vector3.h"
#include "interpolate_pwl.h"


void init_interpolate_pwl(P,Q,p,m)
/* berechnet die Koeffizienten des Interpolationspolynoms */
vector3 	*****P;		/* Koeffizienten des Interpolationspolynoms */
vector3 	***Q;		/* gegebene Knotenpunkte                    */
unsigned int	p; 		/* Anzahl der Pataches  		    */
unsigned int	m;		/* Zahl der Level                           */
{
unsigned int	n = 1<<(m-2);	/* n*n Elemente pro Patch                   */
unsigned int	i1, i2, i3;	/* Laufindizes fuer Elemente                */
unsigned int	j1, j2, j3;	/* Laufindizes fuer Punkte                  */
vector3		q[5], r[25]; 	/* Interpolationspunkte                     */

/* Interpolationsgewichte */
double		a[25] = {     0,     0,     1,     0,      0, \
                          1./12, -2./3,     0,  2./3, -1./12, \
                         -1./24,  2./3, -5./4,  2./3, -1./24, \
                         -1./12,  1./6,     0, -1./6,  1./12, \
			  1./24, -1./6,  1./4, -1./6,  1./24 };

/* Initialisierung */
(*P) = (vector3****) malloc(p*sizeof(vector3***));

for (i1=0; i1<p; i1++)
{  (*P)[i1] = (vector3***) malloc(n*sizeof(vector3**));
   for (i2=0; i2<n; i2++)	/* zeilenweise Durchnumerierung der Patches zi = (i1,i2,i3) */
   {  (*P)[i1][i2] = (vector3**) malloc(n*sizeof(vector3*));
      for (i3=0; i3<n; i3++)
      {  (*P)[i1][i2][i3] = (vector3*) calloc(25,sizeof(vector3));

	 /* Set-up */
         for (j1=0; j1<5; j1++)
	 {  for (j2=0; j2<5; j2++) q[j2] = Q[i1][4*i2+j1][4*i3+j2];
	    for (j2=0; j2<5; j2++)
	    {  r[j2+5*j1].x = r[j2+5*j1].y = r[j2+5*j1].z = 0;
	       for (j3=0; j3<5; j3++)
               {  r[j2+5*j1].x += a[j3+5*j2]*q[j3].x;
                  r[j2+5*j1].y += a[j3+5*j2]*q[j3].y;
                  r[j2+5*j1].z += a[j3+5*j2]*q[j3].z;
                  }
               }
	    }

	 /* berechne Interpolationskoeffizienten */
         for (j1=0; j1<5; j1++)
	 {  for (j2=0; j2<5; j2++)
	    {  for (j3=0; j3<5; j3++)
               {  (*P)[i1][i2][i3][j2+5*j1].x += r[j2+5*j3].x*a[j3+5*j1];
                  (*P)[i1][i2][i3][j2+5*j1].y += r[j2+5*j3].y*a[j3+5*j1];
                  (*P)[i1][i2][i3][j2+5*j1].z += r[j2+5*j3].z*a[j3+5*j1];
                  }
	       }
	    }
         }
      }
   }
return;
}


void free_interpolate_pwl(P,p,m)
/* gibt die Koeffizienten des Interpolationspolynoms frei */
vector3 	*****P;		/* Koeffizienten des Interpolationspolynoms */
unsigned int	p;	 	/* Anzahl der Pataches  		    */
unsigned int	m;		/* Zahl der Level                           */
{
unsigned int	n = 1<<(m-2);	/* n*n Elemente pro Patch                   */
unsigned int	i1, i2, i3;	/* Laufindizes fuer Elemente                */

for (i1=0; i1<p; i1++)
{  for (i2=0; i2<n; i2++)	/* zeilenweise Durchnumerierung der Patches zi = (i1,i2,i3) */
   {  for (i3=0; i3<n; i3++) free((*P)[i1][i2][i3]);
      free((*P)[i1][i2]);
      }
   free((*P)[i1]);
   }
free(*P);
return;
}


vector3 Chi_pwl(a,p,m)
/* wertet das durch p gegebene Interpolationspolynom in a aus */
vector2		a;
vector3		***p;
unsigned int	m;
{
unsigned int	x, y;
vector3		c;

/* Transformation auf [-2,2]^2 */
a.x *= 1<<(m-2);
a.y *= 1<<(m-2);
x = (unsigned int) a.x;
y = (unsigned int) a.y;
a.x = 4*(a.x-x-0.5);
a.y = 4*(a.y-y-0.5);

/* Interpolation */
c.x =       p[y][x][ 0].x+a.x*(p[y][x][ 1].x+a.x*(p[y][x][ 2].x+a.x*(p[y][x][ 3].x+a.x*p[y][x][ 4].x)))    \
    + a.y*((p[y][x][ 5].x+a.x*(p[y][x][ 6].x+a.x*(p[y][x][ 7].x+a.x*(p[y][x][ 8].x+a.x*p[y][x][ 9].x))))   \
    + a.y*((p[y][x][10].x+a.x*(p[y][x][11].x+a.x*(p[y][x][12].x+a.x*(p[y][x][13].x+a.x*p[y][x][14].x))))   \
    + a.y*((p[y][x][15].x+a.x*(p[y][x][16].x+a.x*(p[y][x][17].x+a.x*(p[y][x][18].x+a.x*p[y][x][19].x))))   \
    + a.y* (p[y][x][20].x+a.x*(p[y][x][21].x+a.x*(p[y][x][22].x+a.x*(p[y][x][23].x+a.x*p[y][x][24].x)))))));
c.y =       p[y][x][ 0].y+a.x*(p[y][x][ 1].y+a.x*(p[y][x][ 2].y+a.x*(p[y][x][ 3].y+a.x*p[y][x][ 4].y)))    \
    + a.y*((p[y][x][ 5].y+a.x*(p[y][x][ 6].y+a.x*(p[y][x][ 7].y+a.x*(p[y][x][ 8].y+a.x*p[y][x][ 9].y))))   \
    + a.y*((p[y][x][10].y+a.x*(p[y][x][11].y+a.x*(p[y][x][12].y+a.x*(p[y][x][13].y+a.x*p[y][x][14].y))))   \
    + a.y*((p[y][x][15].y+a.x*(p[y][x][16].y+a.x*(p[y][x][17].y+a.x*(p[y][x][18].y+a.x*p[y][x][19].y))))   \
    + a.y* (p[y][x][20].y+a.x*(p[y][x][21].y+a.x*(p[y][x][22].y+a.x*(p[y][x][23].y+a.x*p[y][x][24].y)))))));
c.z =       p[y][x][ 0].z+a.x*(p[y][x][ 1].z+a.x*(p[y][x][ 2].z+a.x*(p[y][x][ 3].z+a.x*p[y][x][ 4].z)))    \
    + a.y*((p[y][x][ 5].z+a.x*(p[y][x][ 6].z+a.x*(p[y][x][ 7].z+a.x*(p[y][x][ 8].z+a.x*p[y][x][ 9].z))))   \
    + a.y*((p[y][x][10].z+a.x*(p[y][x][11].z+a.x*(p[y][x][12].z+a.x*(p[y][x][13].z+a.x*p[y][x][14].z))))   \
    + a.y*((p[y][x][15].z+a.x*(p[y][x][16].z+a.x*(p[y][x][17].z+a.x*(p[y][x][18].z+a.x*p[y][x][19].z))))   \
    + a.y* (p[y][x][20].z+a.x*(p[y][x][21].z+a.x*(p[y][x][22].z+a.x*(p[y][x][23].z+a.x*p[y][x][24].z)))))));

return(c);
}


vector3 dChi_dx_pwl(a,p,m)
/* definiert die Ableitung nach x des Interpolationspolynoms p */
vector2		a;
vector3		***p;
unsigned int	m;
{
unsigned int	x, y;
vector3		c;
double		w;

/* Transformation auf [-2,2]^2 */
a.x *= 1<<(m-2);
a.y *= 1<<(m-2);
x = (unsigned int) a.x;
y = (unsigned int) a.y;
a.x = 4*(a.x-x-0.5);
a.y = 4*(a.y-y-0.5);
w = 1<<m;

/* Interpolation */
c.x = w*   (p[y][x][ 1].x+a.x*(2*p[y][x][ 2].x+a.x*(3*p[y][x][ 3].x+a.x*4*p[y][x][ 4].x))     \
    + a.y*((p[y][x][ 6].x+a.x*(2*p[y][x][ 7].x+a.x*(3*p[y][x][ 8].x+a.x*4*p[y][x][ 9].x)))    \
    + a.y*((p[y][x][11].x+a.x*(2*p[y][x][12].x+a.x*(3*p[y][x][13].x+a.x*4*p[y][x][14].x)))    \
    + a.y*((p[y][x][16].x+a.x*(2*p[y][x][17].x+a.x*(3*p[y][x][18].x+a.x*4*p[y][x][19].x)))    \
    + a.y* (p[y][x][21].x+a.x*(2*p[y][x][22].x+a.x*(3*p[y][x][23].x+a.x*4*p[y][x][24].x)))))));
c.y = w*   (p[y][x][ 1].y+a.x*(2*p[y][x][ 2].y+a.x*(3*p[y][x][ 3].y+a.x*4*p[y][x][ 4].y))     \
    + a.y*((p[y][x][ 6].y+a.x*(2*p[y][x][ 7].y+a.x*(3*p[y][x][ 8].y+a.x*4*p[y][x][ 9].y)))    \
    + a.y*((p[y][x][11].y+a.x*(2*p[y][x][12].y+a.x*(3*p[y][x][13].y+a.x*4*p[y][x][14].y)))    \
    + a.y*((p[y][x][16].y+a.x*(2*p[y][x][17].y+a.x*(3*p[y][x][18].y+a.x*4*p[y][x][19].y)))    \
    + a.y* (p[y][x][21].y+a.x*(2*p[y][x][22].y+a.x*(3*p[y][x][23].y+a.x*4*p[y][x][24].y)))))));
c.z = w*   (p[y][x][ 1].z+a.x*(2*p[y][x][ 2].z+a.x*(3*p[y][x][ 3].z+a.x*4*p[y][x][ 4].z))     \
    + a.y*((p[y][x][ 6].z+a.x*(2*p[y][x][ 7].z+a.x*(3*p[y][x][ 8].z+a.x*4*p[y][x][ 9].z)))    \
    + a.y*((p[y][x][11].z+a.x*(2*p[y][x][12].z+a.x*(3*p[y][x][13].z+a.x*4*p[y][x][14].z)))    \
    + a.y*((p[y][x][16].z+a.x*(2*p[y][x][17].z+a.x*(3*p[y][x][18].z+a.x*4*p[y][x][19].z)))    \
    + a.y* (p[y][x][21].z+a.x*(2*p[y][x][22].z+a.x*(3*p[y][x][23].z+a.x*4*p[y][x][24].z)))))));

return(c);
}


vector3 dChi_dy_pwl(a,p,m)
/* definiert die Ableitung nach y des Interpolationspolynoms p */
vector2		a;
vector3		***p;
unsigned int	m;
{
unsigned int	x, y;
vector3		c;
double		w;

/* Transformation auf [-2,2]^2 */
a.x *= 1<<(m-2);
a.y *= 1<<(m-2);
x = (unsigned int) a.x;
y = (unsigned int) a.y;
a.x = 4*(a.x-x-0.5);
a.y = 4*(a.y-y-0.5);
w = 1<<m;

/* Interpolation */
c.x = w*  (   p[y][x][ 5].x+a.x*(p[y][x][ 6].x+a.x*(p[y][x][ 7].x+a.x*(p[y][x][ 8].x+a.x*p[y][x][ 9].x)))   \
    + a.y*(2*(p[y][x][10].x+a.x*(p[y][x][11].x+a.x*(p[y][x][12].x+a.x*(p[y][x][13].x+a.x*p[y][x][14].x))))  \
    + a.y*(3*(p[y][x][15].x+a.x*(p[y][x][16].x+a.x*(p[y][x][17].x+a.x*(p[y][x][18].x+a.x*p[y][x][19].x))))  \
    + a.y* 4*(p[y][x][20].x+a.x*(p[y][x][21].x+a.x*(p[y][x][22].x+a.x*(p[y][x][23].x+a.x*p[y][x][24].x)))))));
c.y = w*  (   p[y][x][ 5].y+a.x*(p[y][x][ 6].y+a.x*(p[y][x][ 7].y+a.x*(p[y][x][ 8].y+a.x*p[y][x][ 9].y)))    \
    + a.y*(2*(p[y][x][10].y+a.x*(p[y][x][11].y+a.x*(p[y][x][12].y+a.x*(p[y][x][13].y+a.x*p[y][x][14].y))))   \
    + a.y*(3*(p[y][x][15].y+a.x*(p[y][x][16].y+a.x*(p[y][x][17].y+a.x*(p[y][x][18].y+a.x*p[y][x][19].y))))   \
    + a.y* 4*(p[y][x][20].y+a.x*(p[y][x][21].y+a.x*(p[y][x][22].y+a.x*(p[y][x][23].y+a.x*p[y][x][24].y)))))));
c.z = w*  (   p[y][x][ 5].z+a.x*(p[y][x][ 6].z+a.x*(p[y][x][ 7].z+a.x*(p[y][x][ 8].z+a.x*p[y][x][ 9].z)))    \
    + a.y*(2*(p[y][x][10].z+a.x*(p[y][x][11].z+a.x*(p[y][x][12].z+a.x*(p[y][x][13].z+a.x*p[y][x][14].z))))   \
    + a.y*(3*(p[y][x][15].z+a.x*(p[y][x][16].z+a.x*(p[y][x][17].z+a.x*(p[y][x][18].z+a.x*p[y][x][19].z))))   \
    + a.y* 4*(p[y][x][20].z+a.x*(p[y][x][21].z+a.x*(p[y][x][22].z+a.x*(p[y][x][23].z+a.x*p[y][x][24].z)))))));

return(c);
}


vector3 n_Chi_pwl(a,p,m)
/* definiert die Normalableitung des Interpolationspolynoms p */
vector2		a;
vector3		***p;
unsigned int	m;
{
unsigned int	x, y;
vector3		c, dc_dx, dc_dy;

/* Transformation auf [-2,2]^2 */
a.x *= 1<<(m-2);
a.y *= 1<<(m-2);
x = (unsigned int) a.x;
y = (unsigned int) a.y;
a.x = 4*(a.x-x-0.5);
a.y = 4*(a.y-y-0.5);

/* Interpolation */
dc_dx.x =       p[y][x][ 1].x+a.x*(2*p[y][x][ 2].x+a.x*(3*p[y][x][ 3].x+a.x*4*p[y][x][ 4].x))    \
        + a.y*((p[y][x][ 6].x+a.x*(2*p[y][x][ 7].x+a.x*(3*p[y][x][ 8].x+a.x*4*p[y][x][ 9].x)))   \
        + a.y*((p[y][x][11].x+a.x*(2*p[y][x][12].x+a.x*(3*p[y][x][13].x+a.x*4*p[y][x][14].x)))   \
        + a.y*((p[y][x][16].x+a.x*(2*p[y][x][17].x+a.x*(3*p[y][x][18].x+a.x*4*p[y][x][19].x)))   \
        + a.y* (p[y][x][21].x+a.x*(2*p[y][x][22].x+a.x*(3*p[y][x][23].x+a.x*4*p[y][x][24].x))))));
dc_dx.y =       p[y][x][ 1].y+a.x*(2*p[y][x][ 2].y+a.x*(3*p[y][x][ 3].y+a.x*4*p[y][x][ 4].y))    \
        + a.y*((p[y][x][ 6].y+a.x*(2*p[y][x][ 7].y+a.x*(3*p[y][x][ 8].y+a.x*4*p[y][x][ 9].y)))   \
        + a.y*((p[y][x][11].y+a.x*(2*p[y][x][12].y+a.x*(3*p[y][x][13].y+a.x*4*p[y][x][14].y)))   \
        + a.y*((p[y][x][16].y+a.x*(2*p[y][x][17].y+a.x*(3*p[y][x][18].y+a.x*4*p[y][x][19].y)))   \
        + a.y* (p[y][x][21].y+a.x*(2*p[y][x][22].y+a.x*(3*p[y][x][23].y+a.x*4*p[y][x][24].y))))));
dc_dx.z =       p[y][x][ 1].z+a.x*(2*p[y][x][ 2].z+a.x*(3*p[y][x][ 3].z+a.x*4*p[y][x][ 4].z))    \
        + a.y*((p[y][x][ 6].z+a.x*(2*p[y][x][ 7].z+a.x*(3*p[y][x][ 8].z+a.x*4*p[y][x][ 9].z)))   \
        + a.y*((p[y][x][11].z+a.x*(2*p[y][x][12].z+a.x*(3*p[y][x][13].z+a.x*4*p[y][x][14].z)))   \
        + a.y*((p[y][x][16].z+a.x*(2*p[y][x][17].z+a.x*(3*p[y][x][18].z+a.x*4*p[y][x][19].z)))   \
        + a.y* (p[y][x][21].z+a.x*(2*p[y][x][22].z+a.x*(3*p[y][x][23].z+a.x*4*p[y][x][24].z))))));

dc_dy.x =         p[y][x][ 5].x+a.x*(p[y][x][ 6].x+a.x*(p[y][x][ 7].x+a.x*(p[y][x][ 8].x+a.x*p[y][x][ 9].x)))   \
        + a.y*(2*(p[y][x][10].x+a.x*(p[y][x][11].x+a.x*(p[y][x][12].x+a.x*(p[y][x][13].x+a.x*p[y][x][14].x))))  \
        + a.y*(3*(p[y][x][15].x+a.x*(p[y][x][16].x+a.x*(p[y][x][17].x+a.x*(p[y][x][18].x+a.x*p[y][x][19].x))))  \
        + a.y* 4*(p[y][x][20].x+a.x*(p[y][x][21].x+a.x*(p[y][x][22].x+a.x*(p[y][x][23].x+a.x*p[y][x][24].x))))));
dc_dy.y =         p[y][x][ 5].y+a.x*(p[y][x][ 6].y+a.x*(p[y][x][ 7].y+a.x*(p[y][x][ 8].y+a.x*p[y][x][ 9].y)))   \
        + a.y*(2*(p[y][x][10].y+a.x*(p[y][x][11].y+a.x*(p[y][x][12].y+a.x*(p[y][x][13].y+a.x*p[y][x][14].y))))  \
        + a.y*(3*(p[y][x][15].y+a.x*(p[y][x][16].y+a.x*(p[y][x][17].y+a.x*(p[y][x][18].y+a.x*p[y][x][19].y))))  \
        + a.y* 4*(p[y][x][20].y+a.x*(p[y][x][21].y+a.x*(p[y][x][22].y+a.x*(p[y][x][23].y+a.x*p[y][x][24].y))))));
dc_dy.z =         p[y][x][ 5].z+a.x*(p[y][x][ 6].z+a.x*(p[y][x][ 7].z+a.x*(p[y][x][ 8].z+a.x*p[y][x][ 9].z)))   \
        + a.y*(2*(p[y][x][10].z+a.x*(p[y][x][11].z+a.x*(p[y][x][12].z+a.x*(p[y][x][13].z+a.x*p[y][x][14].z))))  \
        + a.y*(3*(p[y][x][15].z+a.x*(p[y][x][16].z+a.x*(p[y][x][17].z+a.x*(p[y][x][18].z+a.x*p[y][x][19].z))))  \
        + a.y* 4*(p[y][x][20].z+a.x*(p[y][x][21].z+a.x*(p[y][x][22].z+a.x*(p[y][x][23].z+a.x*p[y][x][24].z))))));

c.x = (dc_dx.y*dc_dy.z-dc_dx.z*dc_dy.y)*(1<<2*m);
c.y = (dc_dx.z*dc_dy.x-dc_dx.x*dc_dy.z)*(1<<2*m);
c.z = (dc_dx.x*dc_dy.y-dc_dx.y*dc_dy.x)*(1<<2*m);

return(c);
}
