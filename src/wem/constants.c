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

/*****************
 *  Constants.c  *
 *****************/


/*================================================*
 *  Hier werden alle Paramter gesetzt, die im     *
 *  Wavelet-Galerkin-Verfahren benoetigt werden.  *
 *================================================*/

/* Operatorordnung */
const double op = -1;

/* Kompression: a > 1, 0 < b < 1, d < dp < td-op */

const double a = 2.5;
const double b = 0.001;
const double dp = 2.5;


/* Quadratur */
const unsigned int g_max = 10;  /* maximaler Quadraturgrad             */
const unsigned int min_quadrature_level = 2;    /* minimales Quadraturlevel            */
const double q = 1;             /* Unterteilungskonstante q > 0.25     */
const double scaling_factor = 0.7071;   /* Groesse des relativen Umkreisradius */

/* Genaugikeit der iterativen Loesung */
const double eps = 1e-10;

/* Konstante fuer Feldverlaengerung */
const unsigned int delta = 10;

/* Genaugikeit bei Punktevergleich */
const double tol = 1e-10;
/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

/* warning-disabler-end */

