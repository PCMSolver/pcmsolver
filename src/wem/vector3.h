/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Weffc++"
#pragma clang diagnostic ignored "-Wextra"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdeprecated-register"
#pragma clang diagnostic ignored "-Wincompatible-pointer-types"
#pragma clang diagnostic ignored "-Wempty-body"
#endif

/* warning-disabler-end */

#ifndef VECTOR3
#define VECTOR3
/***************
 *  vector3.h  *
 ***************/


/*====================================================*
 *  Kleine Arithmetik fuer dreidimensionale Vektoren  *
 *====================================================*/


typedef struct {
    double x, y, z;
} vector3;

/* Typdefinition */


vector3 vector3_make(double x, double y, double z);
/* Typkonvertierung: 3xREAL in vector3 */


vector3 vector3_add(vector3 a, vector3 b);
/* Vektoraddition */


vector3 vector3_sub(vector3 a, vector3 b);
/* Vektorsubtraktion */


vector3 vector3_mul(vector3 a, vector3 b);
/* Vektormultiplikation */


vector3 vector3_Smul(double a, vector3 b);
/* S-Multiplikation */


double vector3_skalp(vector3 a, vector3 b);
/* Skalarprodukt */


double vector3_norm(vector3 a);
/* Euklid-Norm */
#endif
/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

/* warning-disabler-end */

