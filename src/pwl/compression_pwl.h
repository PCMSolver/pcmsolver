#ifndef COMPRESSION_pwl
#define COMPRESSION_pwl
/*******************
 *  Compression.h  *
 *******************/


/*===========================================*
 *  Bestimmt in der Steifigkeitsmatrix alle  *
 *  nach der 1. und 2. Kompression zu        *
 *  berechnenden Eintraege.                  *
 *===========================================*/

double compression_pwl(sparse2 *T, wavelet_pwl * W, element_pwl * E, 
                       unsigned int p, unsigned int M, unsigned int np);
#endif
