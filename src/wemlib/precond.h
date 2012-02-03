#ifndef PRECOND
#define PRECOND
/***************
 *  precond.h  *
 ***************/


/*================================================*
 *  Berechnet den Preconditioner und              *
 *  definiert dessen Anwendung auf einen Vektor.  *
 *================================================*/


void inv_G_times_x(double *x, unsigned int p, unsigned int M);
/* Loest das lineare Gleichungssystem T*T'*y = x, wobei x als
   Startvektor verwendet und die Loesung y in x geschrieben wird. */


void precond(double *a, double *b, unsigned int p, unsigned int M);
/* Anwendung des Preconditioners auf den Vektor b, der
   NICHT veraendert wird. Das Ergebnis wird in a gespeichert */
#endif
