#ifndef GAUSS_LEGENDRE
#define GAUSS_LEGENDRE
/**********************
 *  Gauss_Legendre.h  *
 **********************/


/*=======================================*
 *  Stuetzstellen Xi und Gewichte G der  *
 *  Gauss-Quadraturformeln auf [0,1].    *
 *=======================================*/


void            init_Gauss_Legendre(quadrature **Q, unsigned int g);


void            free_Gauss_Legendre(quadrature **Q, unsigned int g);
#endif
