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

#ifndef KERN
#define KERN
/**********
 * kern.h *
 **********/


extern const double epsilon;    /* dielectric constant of the solvent */
extern const double kappa;      /* constant related to ion screening  */

/* inverse dielectric tensor of the solvent */
extern const double epsilon11;
extern const double epsilon12;
extern const double epsilon13;
extern const double epsilon21;
extern const double epsilon22;
extern const double epsilon23;
extern const double epsilon31;
extern const double epsilon32;
extern const double epsilon33;


/*===========================*
 *  Einfachschichtpotential  *
 *===========================*/

double SingleLayerInt(vector3 x, vector3 y);


double SingleLayerExt(vector3 x, vector3 y);


double SingleLayerAni(vector3 x, vector3 y);


/*==========================*
 *  Doppelschichtpotential  *
 *==========================*/

double DoubleLayerInt(vector3 x, vector3 y, vector3 n_y);


double DoubleLayerExt(vector3 x, vector3 y, vector3 n_y);


double DoubleLayerAni(vector3 x, vector3 y, vector3 n_y);
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

