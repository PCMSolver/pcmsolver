#ifndef COMPRESSION
#define COMPRESSION
/*******************
 *  Compression.h  *
 *******************/

/*===========================================*
 *  Bestimmt in der Steifigkeitsmatrix alle  *
 *  nach der 1. und 2. Kompression zu        *
 *  berechnenden Eintraege.                  *
 *===========================================*/


double compression(sparse2 *T, wavelet *W, element *E, unsigned int p, unsigned int M);

#endif
