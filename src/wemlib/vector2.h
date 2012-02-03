#ifndef VECTOR2
#define VECTOR2
/***************
 *  vector2.h  *
 ***************/


/*====================================================*
 *  Kleine Arithmetik fuer zweidimensionale Vektoren  *
 *====================================================*/


typedef struct {
	double          x, y;
} vector2;
/* Typdefinition */


vector2         vector2_make(double x, double y);
/* Typkonvertierung: 2xREAL in vector2 */


vector2         vector2_add(vector2 a, vector2 b);
/* Vektoraddition */


vector2         vector2_sub(vector2 a, vector2 b);
/* Vektorsubtraktion */


vector2         vector2_Smul(double s, vector2 a);
/* S-Multiplikation */


double          vector2_skalp(vector2 a, vector2 b);
/* Skalarprodukt */


double          vector2_norm(vector2 a);
/* Euklid-Norm */
#endif
