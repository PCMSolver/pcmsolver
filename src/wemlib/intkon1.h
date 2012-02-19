#ifndef INTKON1_H_
#define INTKON1_H_
/***************
 *  IntKon1.h  *
 ***************/


/*===========================================================*
 *  Fernfeld-Quadratur-Routine:			             *
 *  Die in der Funktion init_randwerte vorab berechneten     *
 *  Auswertepunkte und Gewichte der Gauss-Quadratur werden   *
 *  in IntKon1 zum entsprechenden Integral zusammengefuegt.  *
 *===========================================================*/


#include "randwerte.h"

void IntKon1(double *c, element *element1, element *element2, randwerte *RW, cubature *Q1, cubature *Q2, vector3 ****P, unsigned int M, double (*SL) (), double (*DL) ());
/* No-Problem-Quadrature-Routine */
#endif
