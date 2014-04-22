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

#ifndef VECTOR2
#define VECTOR2
/***************
 *  vector2.h  *
 ***************/


/*====================================================*
 *  Kleine Arithmetik fuer zweidimensionale Vektoren  *
 *====================================================*/


typedef struct {
    double x, y;
} vector2;

/* Typdefinition */


vector2 vector2_make(double x, double y);
/* Typkonvertierung: 2xREAL in vector2 */


vector2 vector2_add(vector2 a, vector2 b);
/* Vektoraddition */


vector2 vector2_sub(vector2 a, vector2 b);
/* Vektorsubtraktion */


vector2 vector2_Smul(double s, vector2 a);
/* S-Multiplikation */


double vector2_skalp(vector2 a, vector2 b);
/* Skalarprodukt */


double vector2_norm(vector2 a);
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

