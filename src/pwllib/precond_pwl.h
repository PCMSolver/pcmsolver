#ifndef PRECOND_PWL
#define PRECOND_PWL
/***************
 *  precond_pwl.h  *
 ***************/
 

/*================================================*
 *  Berechnet den Preconditioner und              *
 *  definiert dessen Anwendung auf einen Vektor.  *
 *================================================*/
 
 
void inv_A_times_x_pwl(sparse *A, double *x, unsigned int **F, unsigned int p, unsigned int M);
/* Loest das lineare Gleichungssystem T*A*T'*y = x,
   wobei der Nullvektor als Startvektor verwendet 
   und die Loesung y in x geschrieben wird. */


void single_scale_gram_pwl(sparse *G, unsigned int **F, unsigned int p, unsigned int M);
/* berechnet die Massenmatrix in der Einskalenbasis.
   Fuer die Quadratur wird die Mittelpunktsregel verwendet. */


void precond_pwl(double *a, double *b, sparse *G, wavelet_pwl *W, unsigned int **F, unsigned int p, unsigned int M);
/* Anwendung des Preconditioners auf den Vektor b, der
   NICHT veraendert wird. Das Ergebnis wird in a gespeichert */
#endif
