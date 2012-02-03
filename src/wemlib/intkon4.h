#ifndef INTKON4
#define INTKON4
/***************
 *  IntKon4.h  *
 ***************/


void 
IntKon4(double *c, element *element1, element *element2, unsigned int ind_s, unsigned int ind_t,
	cubature *Q, vector3 ****P, unsigned int M, double (*SL) (), double (*DL) ());
/* GEMEINSAME ECKE IM NULLPUNKT */
#endif
