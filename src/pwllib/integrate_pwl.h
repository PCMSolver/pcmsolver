#ifndef INTEGRATE_PWL
#define INTEGRATE_PWL
/*****************
 *  Integrate.h  *
 *****************/
 
 
/*=============================================================*
 *  Enthaelt alle Routinen, die im Wavelet-Galerkin-Verfahren  *
 *  fuer das aufsplitten der Integrale benoetigt werden.       *
 *=============================================================*/
 

void element_element_interaction_pwl(double *c, vector3 *P, element_pwl *E, unsigned int ind1, 
                                     unsigned int ind2, randwerte *RW, cubature *Q, 
                                     vector3 ****R, unsigned int M, double prec, 
                                     double SL(), double DL(), double I);
/* Zerlegungsalgorithmus fuer die Integration Element ind1 mit Element ind2 */
#endif
