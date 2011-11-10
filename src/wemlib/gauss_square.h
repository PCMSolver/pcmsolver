#ifndef GAUSS_SQUARE_H_
#define GAUSS_SQUARE_H_
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
