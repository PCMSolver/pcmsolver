#ifndef QUADRATURE
#define QUADRATURE
/******************
 *  Quadrature.h  *
 ******************/


/*===============================*
 *  Deklariert Quadrtur-Formeln  * 
 *  auf einem Intervall.         *
 *===============================*/


typedef struct 
{  unsigned int	   nop;               /* Zahl der Stuetzstellen und Gewichte         */
   double	   *xi, *w;           /* Stuetzstellen/Gewichte der Quadratur-Formel */
   } quadrature;

#endif
