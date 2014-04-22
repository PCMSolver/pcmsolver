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

#ifndef TRAFOS_PWL
#define TRAFOS_PWL
/**************
 *  Trafos.h  *
 **************/


/*==================================================*
 *  Enthaelt alle Hilfsroutinen fuer die Quadratur  *
 *==================================================*/


vector2 Kappa_pwl(vector2 a, vector2 b, double h);
/* Abbildung vom Einheitsquadrat auf das Element */


vector2 Tau_pwl(double b_x, double b_y, unsigned int CASE);
/* definiert die Drehungen fuer den Duffy-Trick */


unsigned int compare_pwl(vector3 *P, unsigned int *F1, unsigned int *F2, unsigned int *ind1, unsigned int *ind2);
/* Untersucht zwei Patches (F1[0],...,F1[3]) und (F2[0],...,F2[3])
   auf eine gemeinsame Kante (Funktionsergebnis: 3) oder einen 
   gemeinsamen Punkt (Funktionsergebnis: 4), wobei in diesen Faellen 
   die Drehungen ind1 und ind2 entsprechend der Situation richtig gesetzt
   werden. Im Falle keiner Gemeinsamkeiten ist das Funktionsergebnis 1. */


void quadrature_grade_pwl(signed int *g1, signed int *g2, unsigned int m1, unsigned int m2, double dist, double prec);
/* bestimmt fuer Elemente der Level m1 und m2 
   mit dem Abstand dist die benoetigten Quadraturgrade */


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

