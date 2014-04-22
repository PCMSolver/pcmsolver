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
 *  vector3.c  *
 ***************/


/*====================================================*
 *  Kleine Arithmetik fuer dreidimensionale Vektoren  *
 *====================================================*/


#include <math.h>
#include "vector3.h"


vector3 vector3_make(x, y, z)
/* Typkonvertierung: 3xREAL in vector3 */
double x, y, z;
{
    vector3 c;
    c.x = x;
    c.y = y;
    c.z = z;
    return (c);
}


vector3 vector3_add(a, b)
/* Vektoraddition */
vector3 a, b;
{
    vector3 c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return (c);
}


vector3 vector3_sub(a, b)
/* Vektorsubtraktion */
vector3 a, b;
{
    vector3 c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return (c);
}


vector3 vector3_mul(a, b)
/* Vektormultiplikation */
vector3 a, b;
{
    vector3 c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return (c);
}


vector3 vector3_Smul(a, b)
/* S-Multiplikation */
double a;
vector3 b;
{
    vector3 c;
    c.x = a * b.x;
    c.y = a * b.y;
    c.z = a * b.z;
    return (c);
}


double vector3_skalp(a, b)
/* Skalarprodukt */
vector3 a, b;
{
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}


double vector3_norm(a)
/* Euklid-Norm */
vector3 a;
{
    return (sqrt(a.x * a.x + a.y * a.y + a.z * a.z));
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

