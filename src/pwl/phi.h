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

#ifndef PHI
#define PHI
/***********
 *  Phi.h  *
 ***********/


/*============================================*
 *  Definiert die vier stueckweise linearen   *
 *  Basisfunktionen auf dem Einheitsquadrat.  *
 *  Es ist Phi0([0,0]) = 1,		      *
 *         Phi1([1,0]) = 1, 		      *
 * 	   Phi2([1,1]) = 1,		      *
 *         Phi3([0,1]) = 1.		      *
 *============================================*/


/* Ansatzfunktion 0 */
double Phi0(vector2 a);


/* Ansatzfunktion 1 */
double Phi1(vector2 a);


/* Ansatzfunktion 2 */
double Phi2(vector2 a);


/* Ansatzfunktion 3 */
double Phi3(vector2 a);


/* updated c_{i,j} um weight*phi_i(xi)*phi_j(eta) */
void Phi_times_Phi(double *c, double weight, vector2 xi, vector2 eta);


/* updated c_{i,j} um weight * < curl[phi_i(xi)],curl[phi_j(eta)] > */
void Curl_Phi_times_Curl_Phi(double *c, double weight, vector2 xi, vector2 eta, vector3 dChi_dx_pwl_s, vector3 dChi_dy_pwl_s, vector3 dChi_dx_t, vector3 dChi_dy_t);
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

