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

#ifndef DWT_PWL
#define DWT_PWL
/************
 *  dwt1.h  *
 ************/


/*===========================================*
 *  Dieses Modul enthaelt alle Routinen der  *
 *  schnellen Wavelet Transformationen.      *
 *===========================================*/


void multiple(unsigned int ****C, unsigned int **Z, unsigned int **F, unsigned int M, unsigned int p, unsigned int np);
/* Hilfsfunktion zur Wavelet-Transformation: Berechnet aus der 
   Elementliste F eine Liste der lokalen Gitterpunkte (Basisliste) 
   und eine Vielfachheitenliste der Punkte */


void dwtLin(double *a, unsigned int **F, unsigned int M, unsigned int p, unsigned int np);
/* Diskrete Wavelet-Transformation */


void tdwtLin(double *a, unsigned int **F, unsigned int M, unsigned int p, unsigned int np);
/* transformierte Diskrete Wavelet-Transformation */
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

