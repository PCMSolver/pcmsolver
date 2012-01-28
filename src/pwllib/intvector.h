#ifndef INTVECTOR
#define INTVECTOR
/*****************
 *  IntVector.h  *
 *****************/
 

/*================================*
 *  Definiert den Integralvektor  *  
 *================================*/


/*====================*
 *  Typendeklaration  *
 *====================*/
 

typedef struct 
{  
double		sub[48];
}
integral;

typedef struct 
{  
integral	*value;
unsigned int	*index;
unsigned int    integral_number;
}
intvector;


/*===============================================*
 *  Suchalgorithmus gemaess binary-search:       *
 *  Liefert den Index des gewuenschten Elements  *
 *  bzw. -1 falls es nicht vorhanden.            *
 *===============================================*/

signed int search_integral(intvector *I, unsigned int i);


/*=============*
 *  I(i) := z  *
 *=============*/

void set_integral(intvector *I, unsigned int i, double *z);
/* der Eintrag darf nicht vorhanden sein */


/*=======================================*
 *  bestimmt den transponierten Eintrag  *
 *=======================================*/

void permutate(double *a, double *b);
#endif
