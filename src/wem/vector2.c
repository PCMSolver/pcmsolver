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

/***************
 *  vector2.c  *
 ***************/


/*====================================================*
 *  Kleine Arithmetik fuer zweidimensionale Vektoren  *
 *====================================================*/


#include <math.h>
#include "vector2.h"


vector2 vector2_make(x, y)
/* Typkonvertierung: 2xREAL in vector2 */
double x, y;
{
    vector2 c;
    c.x = x;
    c.y = y;
    return (c);
}


vector2 vector2_add(a, b)
/* Vektoraddition */
vector2 a, b;
{
    vector2 c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    return (c);
}


vector2 vector2_sub(a, b)
/* Vektorsubtraktion */
vector2 a, b;
{
    vector2 c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    return (c);
}


vector2 vector2_Smul(a, b)
/* S-Multiplikation */
double a;
vector2 b;
{
    vector2 c;
    c.x = a * b.x;
    c.y = a * b.y;
    return (c);
}


double vector2_skalp(a, b)
/* Skalarprodukt */
vector2 a, b;
{
    return (a.x * b.x + a.y * b.y);
}


double vector2_norm(a)
/* Euklid-Norm */
vector2 a;
{
    return (sqrt(a.x * a.x + a.y * a.y));
}
/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

/* warning-disabler-end */

