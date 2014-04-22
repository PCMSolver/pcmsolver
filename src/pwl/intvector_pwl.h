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

#ifndef INTVECTOR_PWL
#define INTVECTOR_PWL
/*****************
 *  IntVector.h  *
 *****************/


/*================================*
 *  Definiert den Integralvektor  *  
 *================================*/


/*====================*
 *  Typendeklaration  *
 *====================*/


typedef struct {
    double sub[48];
} integral_pwl;

typedef struct {
    integral_pwl *value;
    unsigned int *index;
    unsigned int integral_number;
} intvector_pwl;


/*===============================================*
 *  Suchalgorithmus gemaess binary-search:       *
 *  Liefert den Index des gewuenschten Elements  *
 *  bzw. -1 falls es nicht vorhanden.            *
 *===============================================*/

signed int search_integral_pwl(intvector_pwl * I, unsigned int i);


/*=============*
 *  I(i) := z  *
 *=============*/

void set_integral_pwl(intvector_pwl * I, unsigned int i, double *z);
/* der Eintrag darf nicht vorhanden sein */


/*=======================================*
 *  bestimmt den transponierten Eintrag  *
 *=======================================*/

void permutate(double *a, double *b);
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

