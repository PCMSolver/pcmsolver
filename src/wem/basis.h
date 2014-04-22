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

#ifndef BASIS
#define BASIS
/*************
 *  Basis.h  *
 *************/


/*===========================================*
 *  Hier werden alle Daten bereitgestellt,   *
 *  die fuer das Wavelet-Galerkin-Verfahren  *
 *  benoetigt werden.		             *
 *===========================================*/


typedef struct {                /* Typdefinition Element */
    unsigned int level;         /* Level des Elements                      */
    unsigned int patch;         /* Patch auf dem das Element lebt          */
    unsigned int index_s;       /* Index des Elements in s-Richtung        */
    unsigned int index_t;       /* Index des Elements in t-Richtung        */
    vector3 midpoint;           /* Umkreismittelpunkt                      */
    double radius;              /* Umkreisradius                           */
    unsigned int vertex[4];     /* Indizes der vier Eckpunkte              */
    unsigned int father;        /* Index des Vaterelements                 */
    unsigned int son[4];        /* Indizes der vier Soehne                 */
    unsigned int wavelet_number;        /* Zahl der Wavelets mit diesem Element    */
    unsigned int *wavelet;      /* Indizes der Wavelets mit diesem Element */
    intvector interaction;      /* sparse-Vektor fuer Integrale            */
} element;


typedef struct {                /* Typdefinition Wavelet */
    unsigned int level;         /* Level des Wavelets      */
    unsigned int element_number;        /* number of elements      */
    unsigned int *element;      /* Indizes der Elemente    */
    double *weight;             /* Gewicht der Elemente    */
    unsigned int son[4];        /* Indizes der vier Soehne */
} wavelet;


void unify(vector3 *d, double *r, vector3 d1, double r1, vector3 d2, double r2);
/* bildet die Vereinigung K(d,r) = K(d1,r1) \cup K(d2,r2) */


unsigned int generate_elementlist(element **E, vector3 *P, unsigned int **F, unsigned int p, unsigned int M);
/* erstellt die hierarchsische Elementliste E */


unsigned int generate_waveletlist(wavelet **W, element *E, unsigned int p, unsigned int M);
/* erstellt die Waveletliste W */


void complete_elementlist(wavelet *W, element *E, unsigned int p, unsigned int M);
/* erstellt die Liste zum Zugriff auf die Wavelets */


void simplify_waveletlist(wavelet *W, element *E, unsigned int p, unsigned int M);
/* optimiert die Waveletliste W */


void set_quadrature_level(wavelet *W, element *E, unsigned int p, unsigned int M);
/* verfeinert Grobgitterelemente */


void free_elementlist(element **E, unsigned int p, unsigned int M);
/* gibt den Speicherplatz der hierarchsischen Elementliste E frei */


void free_waveletlist(wavelet **W, unsigned int p, unsigned int M);
/* gibt den Speicherplatz der Waveletliste W frei */


double distance(element *element1, element *element2);
/* Berechnet den Abstand zwischen den Elementen element1 und element2 */
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

