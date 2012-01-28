#ifndef GAUSS_SQUARE
#define GAUSS_SQUARE
/********************
 *  Gauss_Square.h  *
 ********************/
 

/*============================================*
 *  Definiert die Gauss-Quadraturformeln auf  *
 *  dem Referenzviereck [0,1]^2.              *
 *============================================*/
 
 
void init_Gauss_Square(cubature **Q, unsigned int g);


void free_Gauss_Square(cubature **Q, unsigned int g);
#endif
