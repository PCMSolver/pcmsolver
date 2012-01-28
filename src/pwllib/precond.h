#ifndef PRECOND
#define PRECOND
/***************
 *  precond.h  *
 ***************/
 

/*================================================*
 *  Berechnet den Preconditioner und              *
 *  definiert dessen Anwendung auf einen Vektor.  *
 *================================================*/
 
 
void inv_A_times_x(sparse *A, double *x, unsigned int **F, unsigned int p, unsigned int M);
/* Loest das lineare Gleichungssystem T*A*T'*y = x,
   wobei der Nullvektor als Startvektor verwendet 
   und die Loesung y in x geschrieben wird. */


void single_scale_gram(sparse *G, unsigned int **F, unsigned int p, unsigned int M);
/* berechnet die Massenmatrix in der Einskalenbasis.
   Fuer die Quadratur wird die Mittelpunktsregel verwendet. */


void precond(double *a, double *b, sparse *G, wavelet *W, unsigned int **F, unsigned int p, unsigned int M);
/* Anwendung des Preconditioners auf den Vektor b, der
   NICHT veraendert wird. Das Ergebnis wird in a gespeichert */
#endif
