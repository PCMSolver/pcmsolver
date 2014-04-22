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

#ifndef WEMPCG_H_
#define WEMPCG_H_
/**************
 *  WEMPCG.h  *
 **************/


/*=============================================================*
 *  WEMPCG(A,b,x,epsi,p,M)                                     *
 *	                                                       *
 *  Verfahren der konjugierten Gradienten zur Loesung des      *
 *  linearen Gleichungssystems                                 *
 *                                                             *
 *            A1*x = -(G*A2^(-1)/(epsilon-1)+I)*b.             *
 *                                                             *
 *  Vorkonditionierung per Wavelet-Preconditioner.             *
 *	                                                       *
 *  Parameter :                                                *
 *		A    : Matrix im sparse2-Format                *
 *		G    : Massenmatrix                            *
 *		b    : rechte Seite                            *
 *		x    : Startwert und Endwert                   *
 *		epsi : Genauigkeit                             *
 *		p    : Anzahl der Paramtergebiete              *
 *		M    : Zahl der Level                          *
 *=============================================================*/


unsigned int WEMPCG(sparse2 *A, double *b, double *x, double epsi, unsigned int p, unsigned int M);
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

