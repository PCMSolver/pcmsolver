/************
 *  kern.c  *
 ************/


#include <math.h>
#include "vector2.h"
#include "vector3.h"
#include "kern.h"


const double	epsilon = 78.39; 	/* dielectric constant of the solvent */ 
const double	kappa = 0.0;		/* constant related to ion screening  */

/* inverse dielectric tensor of the solvent */ 
const double	epsilon11 = 0.01275672917464;
const double	epsilon12 = 0;
const double	epsilon13 = 0;
const double	epsilon21 = 0;
const double	epsilon22 = 0.01275672917464;
const double	epsilon23 = 0;
const double	epsilon31 = 0;
const double	epsilon32 = 0;
const double	epsilon33 = 0.01275672917464;

/*===========================*
 *  Einfachschichtpotential  *
 *===========================*/

double SingleLayerInt(x,y)
vector3		x, y;
{  
return(1/sqrt((x.x-y.x)*(x.x-y.x)+(x.y-y.y)*(x.y-y.y)+(x.z-y.z)*(x.z-y.z)));
}


double SingleLayerExt(x,y)
vector3		x, y;
{
double 		r = sqrt((x.x-y.x)*(x.x-y.x)+(x.y-y.y)*(x.y-y.y)+(x.z-y.z)*(x.z-y.z));
return(exp(-kappa*r)/(r*epsilon));
}


double SingleLayerAni(x,y)
vector3		x, y;
{
vector3		c;
double		r, det;
c.x = x.x-y.x;
c.y = x.y-y.y;
c.z = x.z-y.z;
det = sqrt( epsilon11*epsilon22*epsilon33 + epsilon12*epsilon23*epsilon31 + epsilon13*epsilon21*epsilon32 \
          - epsilon11*epsilon32*epsilon23 - epsilon21*epsilon12*epsilon33 - epsilon31*epsilon22*epsilon13 );
r = sqrt( epsilon11*c.x*c.x+epsilon12*c.x*c.y+epsilon13*c.x*c.z \
        + epsilon21*c.y*c.x+epsilon22*c.y*c.y+epsilon23*c.y*c.z \
        + epsilon31*c.z*c.x+epsilon32*c.z*c.y+epsilon33*c.z*c.z );
return(det/r);
}


/*==========================*
 *  Doppelschichtpotential  *
 *==========================*/

double DoubleLayerInt(x,y,n_y)
vector3		x, y, n_y;
{  
vector3		c;
double		r;
c.x = x.x-y.x;
c.y = x.y-y.y;
c.z = x.z-y.z;
r = sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
return((c.x*n_y.x+c.y*n_y.y+c.z*n_y.z)/(r*r*r));
}


double DoubleLayerExt(x,y,n_y)
vector3		x, y, n_y;
{  
vector3		c;
double		r;
c.x = x.x-y.x;
c.y = x.y-y.y;
c.z = x.z-y.z;
r = sqrt(c.x*c.x+c.y*c.y+c.z*c.z);
return(exp(-kappa*r)*(c.x*n_y.x+c.y*n_y.y+c.z*n_y.z)*(1+kappa*r)/(r*r*r));
/*
return(exp(-kappa*r)*(c.x*n_y.x+c.y*n_y.y+c.z*n_y.z)*(1+kappa*r)/(r*r*r)/epsilon);
/* /epsilon ? */
}


double DoubleLayerAni(x,y,n_y)
vector3		x, y, n_y;
{  
vector3		c;
double		r, det;
c.x = x.x-y.x;
c.y = x.y-y.y;
c.z = x.z-y.z;
det = sqrt( epsilon11*epsilon22*epsilon33 + epsilon12*epsilon23*epsilon31 + epsilon13*epsilon21*epsilon32 \
          - epsilon11*epsilon32*epsilon23 - epsilon21*epsilon12*epsilon33 - epsilon31*epsilon22*epsilon13 );
r = sqrt( epsilon11*c.x*c.x+epsilon12*c.x*c.y+epsilon13*c.x*c.z \
        + epsilon21*c.y*c.x+epsilon22*c.y*c.y+epsilon23*c.y*c.z \
        + epsilon31*c.z*c.x+epsilon32*c.z*c.y+epsilon33*c.z*c.z );
return(det*(c.x*n_y.x+c.y*n_y.y+c.z*n_y.z)/(r*r*r));
}
