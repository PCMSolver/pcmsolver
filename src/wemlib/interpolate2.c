/*******************
 *  Interpolate.c  *
 *******************/


/*=============================================================*
 *  Enthaelt alle Routinen zur Interpolation der Oberflaeche.  *
 *=============================================================*/
 

#include <stdlib.h>
#include "vector2.h"
#include "vector3.h"
#include "interpolate.h"


void init_interpolate(P,Q,p,m)
/* berechnet die Koeffizienten des Interpolationspolynoms */
vector3 	*****P;		/* Koeffizienten des Interpolationspolynoms */
vector3 	***Q;		/* gegebene Knotenpunkte                    */
unsigned int	p; 		/* Anzahl der Pataches  		    */
unsigned int	m;		/* Zahl der Level                           */
{
unsigned int	n = 1<<(m-1);	/* n*n Elemente pro Patch                   */
unsigned int	i1, i2, i3;	/* Laufindizes fuer Elemente                */
vector3		q[9]; 		/* Interpolationspunkte                     */

/* Initialisierung */
(*P) = (vector3****) malloc(p*sizeof(vector3***));

for (i1=0; i1<p; i1++)
{  (*P)[i1] = (vector3***) malloc(n*sizeof(vector3**));
   for (i2=0; i2<n; i2++)	/* zeilenweise Durchnumerierung der Patches zi = (i1,i2,i3) */
   {  (*P)[i1][i2] = (vector3**) malloc(n*sizeof(vector3*));
      for (i3=0; i3<n; i3++) 	
      {  (*P)[i1][i2][i3] = (vector3*) malloc(9*sizeof(vector3));
	 
	 /* bestimme Interpolationspunkte */
	 q[0] = Q[i1][2*i2  ][2*i3  ];
	 q[1] = Q[i1][2*i2  ][2*i3+1];
	 q[2] = Q[i1][2*i2  ][2*i3+2];
	 q[3] = Q[i1][2*i2+1][2*i3  ];
	 q[4] = Q[i1][2*i2+1][2*i3+1];
	 q[5] = Q[i1][2*i2+1][2*i3+2];
	 q[6] = Q[i1][2*i2+2][2*i3  ];
	 q[7] = Q[i1][2*i2+2][2*i3+1];
	 q[8] = Q[i1][2*i2+2][2*i3+2];
	 
	 /* berechne Koeffizienten */
	 (*P)[i1][i2][i3][0].x = 4*(q[0].x-2*q[1].x+q[2].x-2*q[3].x+4*q[4].x-2*q[5].x+q[6].x-2*q[7].x+q[8].x);
	 (*P)[i1][i2][i3][0].y = 4*(q[0].y-2*q[1].y+q[2].y-2*q[3].y+4*q[4].y-2*q[5].y+q[6].y-2*q[7].y+q[8].y);
	 (*P)[i1][i2][i3][0].z = 4*(q[0].z-2*q[1].z+q[2].z-2*q[3].z+4*q[4].z-2*q[5].z+q[6].z-2*q[7].z+q[8].z);

	 (*P)[i1][i2][i3][1].x = 2*(-3*q[0].x+4*q[1].x-q[2].x+6*q[3].x-8*q[4].x+2*q[5].x-3*q[6].x+4*q[7].x-q[8].x);
	 (*P)[i1][i2][i3][1].y = 2*(-3*q[0].y+4*q[1].y-q[2].y+6*q[3].y-8*q[4].y+2*q[5].y-3*q[6].y+4*q[7].y-q[8].y);
	 (*P)[i1][i2][i3][1].z = 2*(-3*q[0].z+4*q[1].z-q[2].z+6*q[3].z-8*q[4].z+2*q[5].z-3*q[6].z+4*q[7].z-q[8].z);

	 (*P)[i1][i2][i3][2].x = 2*(q[0].x+q[6].x-2*q[3].x);
	 (*P)[i1][i2][i3][2].y = 2*(q[0].y+q[6].y-2*q[3].y);
	 (*P)[i1][i2][i3][2].z = 2*(q[0].z+q[6].z-2*q[3].z);

 	 (*P)[i1][i2][i3][3].x = 2*(-3*q[0].x+6*q[1].x-3*q[2].x+4*q[3].x-8*q[4].x+4*q[5].x-q[6].x+2*q[7].x-q[8].x);
	 (*P)[i1][i2][i3][3].y = 2*(-3*q[0].y+6*q[1].y-3*q[2].y+4*q[3].y-8*q[4].y+4*q[5].y-q[6].y+2*q[7].y-q[8].y);
	 (*P)[i1][i2][i3][3].z = 2*(-3*q[0].z+6*q[1].z-3*q[2].z+4*q[3].z-8*q[4].z+4*q[5].z-q[6].z+2*q[7].z-q[8].z);

	 (*P)[i1][i2][i3][4].x = 9*q[0].x-12*q[1].x+3*q[2].x-12*q[3].x+16*q[4].x-4*q[5].x+3*q[6].x-4*q[7].x+q[8].x;
	 (*P)[i1][i2][i3][4].y = 9*q[0].y-12*q[1].y+3*q[2].y-12*q[3].y+16*q[4].y-4*q[5].y+3*q[6].y-4*q[7].y+q[8].y;
	 (*P)[i1][i2][i3][4].z = 9*q[0].z-12*q[1].z+3*q[2].z-12*q[3].z+16*q[4].z-4*q[5].z+3*q[6].z-4*q[7].z+q[8].z;

	 (*P)[i1][i2][i3][5].x = 4*q[3].x-3*q[0].x-q[6].x;
	 (*P)[i1][i2][i3][5].y = 4*q[3].y-3*q[0].y-q[6].y;
	 (*P)[i1][i2][i3][5].z = 4*q[3].z-3*q[0].z-q[6].z;

	 (*P)[i1][i2][i3][6].x = 2*(q[0].x+q[2].x-2*q[1].x);
	 (*P)[i1][i2][i3][6].y = 2*(q[0].y+q[2].y-2*q[1].y);
	 (*P)[i1][i2][i3][6].z = 2*(q[0].z+q[2].z-2*q[1].z);

	 (*P)[i1][i2][i3][7].x = 4*q[1].x-3*q[0].x-q[2].x;
	 (*P)[i1][i2][i3][7].y = 4*q[1].y-3*q[0].y-q[2].y;
	 (*P)[i1][i2][i3][7].z = 4*q[1].z-3*q[0].z-q[2].z;

	 (*P)[i1][i2][i3][8].x = q[0].x;
	 (*P)[i1][i2][i3][8].y = q[0].y;
	 (*P)[i1][i2][i3][8].z = q[0].z;
	 }
      }
   }
return;
}


void free_interpolate(P,p,m)
/* gibt die Koeffizienten des Interpolationspolynoms frei */
vector3 	*****P;		/* Koeffizienten des Interpolationspolynoms */
unsigned int	p; 		/* Anzahl der Pataches  		    */
unsigned int	m;		/* Zahl der Level                           */
{
unsigned int	n = 1<<(m-1);	/* n*n Elemente pro Patch                   */
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


vector3 Chi(a,p,m)
/* wertet das durch p gegebene Interpolationspolynom in a aus */
vector2		a;
vector3		***p;
unsigned int	m;
{
unsigned int	x, y;
vector3		c;

/* Transfomration auf [0,1]^2 */
a.x *= 1<<(m-1);
a.y *= 1<<(m-1);
x = (unsigned int) a.x;
y = (unsigned int) a.y;
a.x -= x;
a.y -= y;

/* Interpolation */
c.x = (a.y*((p[y][x][0].x*a.x + p[y][x][1].x)*a.x + p[y][x][2].x) + (p[y][x][3].x*a.x + p[y][x][4].x)*a.x + p[y][x][5].x)*a.y + (p[y][x][6].x*a.x + p[y][x][7].x)*a.x + p[y][x][8].x;
c.y = (a.y*((p[y][x][0].y*a.x + p[y][x][1].y)*a.x + p[y][x][2].y) + (p[y][x][3].y*a.x + p[y][x][4].y)*a.x + p[y][x][5].y)*a.y + (p[y][x][6].y*a.x + p[y][x][7].y)*a.x + p[y][x][8].y;
c.z = (a.y*((p[y][x][0].z*a.x + p[y][x][1].z)*a.x + p[y][x][2].z) + (p[y][x][3].z*a.x + p[y][x][4].z)*a.x + p[y][x][5].z)*a.y + (p[y][x][6].z*a.x + p[y][x][7].z)*a.x + p[y][x][8].z;

return(c);
}


vector3 dChi_dx(a,p,m)
/* definiert die Ableitung nach x des Interpolationspolynoms p */
vector2		a;
vector3		***p;
unsigned int	m;
{
unsigned int	x, y;
vector3		c;

/* Transfomration auf [0,1]^2 */
a.x *= 1<<(m-1);
a.y *= 1<<(m-1);
x = (unsigned int) a.x;
y = (unsigned int) a.y;
a.x -= x;
a.y -= y;

/* Interpolation */
c.x = (a.y*(a.y*(2*p[y][x][0].x*a.x + p[y][x][1].x) + 2*p[y][x][3].x*a.x + p[y][x][4].x) + 2*p[y][x][6].x*a.x + p[y][x][7].x)*(1<<(m-1));
c.y = (a.y*(a.y*(2*p[y][x][0].y*a.x + p[y][x][1].y) + 2*p[y][x][3].y*a.x + p[y][x][4].y) + 2*p[y][x][6].y*a.x + p[y][x][7].y)*(1<<(m-1));
c.z = (a.y*(a.y*(2*p[y][x][0].z*a.x + p[y][x][1].z) + 2*p[y][x][3].z*a.x + p[y][x][4].z) + 2*p[y][x][6].z*a.x + p[y][x][7].z)*(1<<(m-1));

return(c);
}


vector3 dChi_dy(a,p,m)
/* definiert die Ableitung nach y des Interpolationspolynoms p */
vector2		a;
vector3		***p;
unsigned int	m;
{
unsigned int	x, y;
vector3		c;

/* Transfomration auf [0,1]^2 */
a.x *= 1<<(m-1);
a.y *= 1<<(m-1);
x = (unsigned int) a.x;
y = (unsigned int) a.y;
a.x -= x;
a.y -= y;

/* Interpolation */
c.x = (2*a.y*((p[y][x][0].x*a.x + p[y][x][1].x)*a.x + p[y][x][2].x) + (p[y][x][3].x*a.x + p[y][x][4].x)*a.x + p[y][x][5].x)*(1<<(m-1));
c.y = (2*a.y*((p[y][x][0].y*a.x + p[y][x][1].y)*a.x + p[y][x][2].y) + (p[y][x][3].y*a.x + p[y][x][4].y)*a.x + p[y][x][5].y)*(1<<(m-1));
c.z = (2*a.y*((p[y][x][0].z*a.x + p[y][x][1].z)*a.x + p[y][x][2].z) + (p[y][x][3].z*a.x + p[y][x][4].z)*a.x + p[y][x][5].z)*(1<<(m-1));

return(c);
}


vector3 n_Chi(a,p,m)
/* definiert die Normalableitung des Interpolationspolynoms p */
vector2		a;
vector3		***p;
unsigned int	m;
{
unsigned int	x, y;
vector3		c, dc_dx, dc_dy;

/* Transfomration auf [0,1]^2 */
a.x *= 1<<(m-1);
a.y *= 1<<(m-1);
x = (unsigned int) a.x;
y = (unsigned int) a.y;
a.x -= x;
a.y -= y;

/* Interpolation */
dc_dx.x = a.y*(a.y*(2*p[y][x][0].x*a.x + p[y][x][1].x) + 2*p[y][x][3].x*a.x + p[y][x][4].x) + 2*p[y][x][6].x*a.x + p[y][x][7].x;
dc_dx.y = a.y*(a.y*(2*p[y][x][0].y*a.x + p[y][x][1].y) + 2*p[y][x][3].y*a.x + p[y][x][4].y) + 2*p[y][x][6].y*a.x + p[y][x][7].y;
dc_dx.z = a.y*(a.y*(2*p[y][x][0].z*a.x + p[y][x][1].z) + 2*p[y][x][3].z*a.x + p[y][x][4].z) + 2*p[y][x][6].z*a.x + p[y][x][7].z;
dc_dy.x = 2*a.y*((p[y][x][0].x*a.x + p[y][x][1].x)*a.x + p[y][x][2].x) + (p[y][x][3].x*a.x + p[y][x][4].x)*a.x + p[y][x][5].x;
dc_dy.y = 2*a.y*((p[y][x][0].y*a.x + p[y][x][1].y)*a.x + p[y][x][2].y) + (p[y][x][3].y*a.x + p[y][x][4].y)*a.x + p[y][x][5].y;
dc_dy.z = 2*a.y*((p[y][x][0].z*a.x + p[y][x][1].z)*a.x + p[y][x][2].z) + (p[y][x][3].z*a.x + p[y][x][4].z)*a.x + p[y][x][5].z;
c.x = (dc_dx.y*dc_dy.z - dc_dx.z*dc_dy.y)*(1<<2*(m-1));
c.y = (dc_dx.z*dc_dy.x - dc_dx.x*dc_dy.z)*(1<<2*(m-1));
c.z = (dc_dx.x*dc_dy.y - dc_dx.y*dc_dy.x)*(1<<2*(m-1));

return(c);
}
