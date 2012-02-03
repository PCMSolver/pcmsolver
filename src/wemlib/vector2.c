/***************
 *  vector2.c  *
 ***************/


/*====================================================*
 *  Kleine Arithmetik fuer zweidimensionale Vektoren  *
 *====================================================*/


#include <math.h>
#include "vector2.h"


vector2 
vector2_make(x, y)
/* Typkonvertierung: 2xREAL in vector2 */
	double          x, y;
{
	vector2         c;
	c.x = x;
	c.y = y;
	return (c);
}


vector2 
vector2_add(a, b)
/* Vektoraddition */
	vector2         a, b;
{
	vector2         c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	return (c);
}


vector2 
vector2_sub(a, b)
/* Vektorsubtraktion */
	vector2         a, b;
{
	vector2         c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	return (c);
}


vector2 
vector2_Smul(a, b)
/* S-Multiplikation */
	double          a;
	vector2         b;
{
	vector2         c;
	c.x = a * b.x;
	c.y = a * b.y;
	return (c);
}


double 
vector2_skalp(a, b)
/* Skalarprodukt */
	vector2         a, b;
{
	return (a.x * b.x + a.y * b.y);
}


double 
vector2_norm(a)
/* Euklid-Norm */
	vector2         a;
{
	return (sqrt(a.x * a.x + a.y * a.y));
}
