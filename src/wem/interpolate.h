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

#ifndef INTERPOLATE
#define INTERPOLATE
/*******************
 *  Interpolate.h  *
 *******************/
#ifdef __cplusplus
extern "C" {
#endif

#include "vector2.h"

/*=============================================================*
 *  Enthaelt alle Routinen zur Interpolation der Oberflaeche.  *
 *=============================================================*/


    void init_interpolate(vector3 *****P, vector3 ***Q, unsigned int p, unsigned int m);
/* berechnet die Koeffizienten des Interpolationspolynoms */


    void free_interpolate(vector3 *****P, unsigned int p, unsigned int m);
/* gibt die Koeffizienten des Interpolationspolynoms frei */


    vector3 Chi(vector2 a, vector3 ***p, unsigned int m);
/* wertet das durch p gegebene Interpolationspolynom in a aus */


    vector3 dChi_dx(vector2 a, vector3 ***p, unsigned int m);
/* definiert die Ableitung nach x des Interpolationspolynoms p */


    vector3 dChi_dy(vector2 a, vector3 ***p, unsigned int m);
/* definiert die Ableitung nach y des Interpolationspolynoms p */


    vector3 n_Chi(vector2 a, vector3 ***p, unsigned int m);
/* definiert die Normalableitung des Interpolationspolynoms p */

#ifdef __cplusplus
}
#endif
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

