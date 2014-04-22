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

#ifndef CONSTANTS
#define CONSTANTS
/*****************
 *  Constants.h  *
 *****************/


/*================================================*
 *  Hier werden alle Paramter gesetzt, die im     *
 *  Wavelet-Galerkin-Verfahren benoetigt werden.  *
 *================================================*/

/* Operatorordnung */
extern const double op;

/* Konstanten bezueglich der verwendeteten Wavelet-Basis */
extern const unsigned int td;
extern const unsigned int minLevel;
extern const unsigned int td_pwl;
extern const unsigned int minLevel_pwl;

/* Kompression: a > 1, 0 < b < 1, d < dp < td-op */
extern const double a;
extern const double b;
extern const double dp;

/* Quadratur */
extern const unsigned int g_max;        /* maximaler Quadraturgrad             */
extern const unsigned int min_quadrature_level; /* minimales Quadraturlevel            */
extern const double q;          /* Unterteilungskonstante q > 0.25     */
extern const double scaling_factor;     /* Groesse des relativen Umkreisradius */

/* Genaugikeit der iterativen Loesung */
extern const double eps;

/* Konstante fuer Feldverlaengerung */
extern const unsigned int delta;

/* Genaugikeit bei Punktevergleich */
extern const double tol;
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

