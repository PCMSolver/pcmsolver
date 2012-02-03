#ifndef VECTOR3
#define VECTOR3
/***************
 *  vector3.h  *
 ***************/


/*====================================================*
 *  Kleine Arithmetik fuer dreidimensionale Vektoren  *
 *====================================================*/


typedef struct {
	double          x, y, z;
} vector3;
/* Typdefinition */


vector3         vector3_make(double x, double y, double z);
/* Typkonvertierung: 3xREAL in vector3 */


vector3         vector3_add(vector3 a, vector3 b);
/* Vektoraddition */


vector3         vector3_sub(vector3 a, vector3 b);
/* Vektorsubtraktion */


vector3         vector3_mul(vector3 a, vector3 b);
/* Vektormultiplikation */


vector3         vector3_Smul(double a, vector3 b);
/* S-Multiplikation */


double          vector3_skalp(vector3 a, vector3 b);
/* Skalarprodukt */


double          vector3_norm(vector3 a);
/* Euklid-Norm */
#endif
