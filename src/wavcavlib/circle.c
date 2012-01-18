/*
 * Purpose: Patch representation of molecular cavities from
 *			atomic coordinates and radii.
 * Version: September 21, 2009.
 * Author : Maharavo Randrianarivony.
 * Affiliation: Institut für Numerische Simulation.
 *              Rheinische Friedrich-Wilhelm Universität Bonn.
 *              Wegelerstraße 6, Bonn 53115.
 *              Germany.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"


int nujw_seco_lijm(double a,double b,double c,double *x)
{int sol;
double dl,rac,eps=1.0e-10,den;
dl=b*b-4*a*c;
if(dl<0)
	sol=0;
else 
	{rac=sqrt(dl);
	den=2.0*a;
	x[0]=(-b+rac)/den;
	x[1]=(-b-rac)/den;
	if(rac<eps)
		sol=1;
	else
		sol=2;
	}
return sol;
}



int qudl_circ_vufw(parm omega,double r,
double a,double b,double c,parm *res)
{int sol,i;
double x_om,y_om,*x,A,B,C,B1,B2,B3,C1,C2,C3,den;
x_om=omega.u;
y_om=omega.v;
x=(double *)malloc(2*sizeof(double));
if(fabs(a)<fabs(b))
	{den=b*b;
	A=1.0+((a*a)/den);
	B1=2*a*c/den;
	B2=2*a*y_om/b;
	B3=-2*x_om;
	B=B1+B2+B3;
	C1=(c*c)/den;
	C2=x_om*x_om+y_om*y_om;
	C3=2*y_om*c/b;
	C=C1+C2+C3-(r*r);
	sol=nujw_seco_lijm(A,B,C,x);
	if(sol!=0)
		{for(i=0;i<2;i++)
			{res[i].u=x[i];
			res[i].v=(-a*x[i]-c)/b;
			}
		}
	}
else
	{den=a*a;
	A=1.0+((b*b)/den);
	B1=2*b*c/den;
	B2=2*b*x_om/a;
	B3=-2*y_om;
	B=B1+B2+B3;
	C1=(c*c)/den;
	C2=x_om*x_om+y_om*y_om;
	C3=2*x_om*c/a;
	C=C1+C2+C3-(r*r);
	sol=nujw_seco_lijm(A,B,C,x);
	if(sol!=0)
		{for(i=0;i<2;i++)
			{res[i].u=(-b*x[i]-c)/a;
			res[i].v=x[i];
			}
		}
	}
free(x);
return sol;
}



int sulw_circ_mojt(parm omega1,double r1,
parm omega2,double r2,parm *res)
{int sol;
double x1,x2,y1,y2,A,B,C;
x1=omega1.u;	y1=omega1.v;
x2=omega2.u;	y2=omega2.v;
A=-2.0*(x1-x2);
B=-2.0*(y1-y2);
C=(y1*y1-y2*y2)+(x1*x1-x2*x2)+(r2*r2-r1*r1);
sol=qudl_circ_vufw(omega1,r1,A,B,C,res);
return sol;
}




