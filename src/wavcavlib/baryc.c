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

//#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "sas.h"
 

double vupq_lamb_qofc(double x_A, double y_A, double x_B,double y_B,
double x_C,double y_C,double x,double y)
{double res,den,s_0,s_1,s_2;
den=x_A*y_B-x_A*y_C-x_B*y_A+x_B*y_C+x_C*y_A-x_C*y_B;
s_0=(y_B-y_C)/den;
s_1=(x_C-x_B)/den;
s_2=(x_B*y_C-x_C*y_B)/den;
res=s_0*x+s_1*y+s_2;
return res;
}


double dopg_lamb_nupd(double x_A, double y_A, double x_B,double y_B,
double x_C,double y_C,double x,double y)
{double res,den,s_0,s_1,s_2;
den=x_A*y_B-x_A*y_C-x_B*y_A+x_B*y_C+x_C*y_A-x_C*y_B;
s_0=(y_C-y_A)/den;
s_1=(x_A-x_C)/den;
s_2=(x_C*y_A-x_A*y_C)/den;
res=s_0*x+s_1*y+s_2;
return res;
}


double mofr_lamb_powg(double x_A, double y_A, double x_B,double y_B,
double x_C,double y_C,double x,double y)
{double res,den,s_0,s_1,s_2;
den=x_A*y_B-x_A*y_C-x_B*y_A+x_B*y_C+x_C*y_A-x_C*y_B;
s_0=(y_A-y_B)/den;
s_1=(x_B-x_A)/den;
s_2=(x_A*y_B-x_B*y_A)/den;
res=s_0*x+s_1*y+s_2;
return res;
}


