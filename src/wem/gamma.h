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

#ifndef GAMMA
#define GAMMA
/*************
 *  Gamma.h  *
 *************/


/*===========================================*
 *  Dummy-Header fuer die Parametrisierung,  *
 *  im makefile ist dann das entsprechende   *
 *  c-File zu setzen.			     *
 *===========================================*/


typedef struct {
    vector3 (*f) (vector2 a);   /* Parametrisierung     */
    vector3 (*df_dx) (vector2 a);       /* Ableitung nach x     */
    vector3 (*df_dy) (vector2 a);       /* Ableitung nach y     */
    vector3 (*n_f) (vector2 a); /* zugehoerende Normale */
} parametrix;


unsigned int init_p(void);
/* Initialisierung: Liefert als Funktionsergebnis die
   Anzahl p der Parametergebiete */


void init_Chi(parametrix **Chi);
/* allokiert den noetigen Speicherplatz fuer die
   Parametrisierung und definiert Chi[0],...,Chi[p] */


void free_Chi(parametrix **Chi);
/* gibt den Speicherplatz fuer die Parametrisierung frei */
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

