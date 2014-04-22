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

#ifndef WEMPGMRES_pwl
#define WEMPGMRES_pwl
/*****************
 *  WEMPGMRES_pwl.h  *
 *****************/


/*==============================================================*
 *  WEMPGMRES_pwl(A,b,x,epsi,W,F,p,M)                               *
 *	                                                        *
 *  GMRES-Verfahren zur Loesung des linearen Gleichungssystems  *
 *	                                                        *
 *		 A2'*x = b  bzw. A2*x = b.                      *
 *	                                                        *
 *  Vorkonditionierung per Diagonalskalierung.                  *
 *	                                                        *
 *  Parameter :                                                 *
 *		A    : Matrix im sparse2-Format                 *
 *		b    : rechte Seite                             *
 *		x    : Startwert und Endwert                    *
 *		epsi : Genauigkeit                              *
 *		W    : Liste der Wavelets                       *
 *		F    : Elementliste der Einskalenbasis          *
 *		p    : Anzahl der Paramtergebiete               *
 *		M    : Zahl der Level                           *
 *==============================================================*/


unsigned int WEMPGMRES1_pwl(sparse2 *A, double *b, double *x, double epsi, wavelet_pwl * W, unsigned int **F, unsigned int p, unsigned int M);


unsigned int WEMPGMRES2_pwl(sparse2 *A, double *b, double *x, double epsi, wavelet_pwl * W, unsigned int **F, unsigned int p, unsigned int M);


/*==============================================================*
 *  WEMPGMRES_pwl(A,B,rhs,x,epsi,W,F,p,M)                           *
 *	                                                        *
 *  GMRES-Verfahren zur Loesung des linearen Gleichungssystems  *
 *	                                                        *
 *	    (B1*G^(-1)*A2'-B2*G^(-1)*A1)*x = rhs.               *
 *	                                                        *
 *  Vorkonditionierung per Wavelet-Preconditioner.              *
 *	                                                        *
 *  Parameter :                                                 *
 *		A, B : Matrizen im sparse2-Format               *
 *		rhs  : rechte Seite                             *
 *		x    : Startwert und Endwert                    *
 *		epsi : Genauigkeit                              *
 *		W    : Liste der Wavelets                       *
 *		F    : Elementliste der Einskalenbasis          *
 *		p    : Anzahl der Paramtergebiete               *
 *		M    : Zahl der Level                           *
 *==============================================================*/

unsigned int WEMPGMRES3_pwl(sparse2 *A, sparse2 *B, double *rhs, double *x, double epsi, wavelet_pwl * W, unsigned int **F, unsigned int p, unsigned int M);
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

