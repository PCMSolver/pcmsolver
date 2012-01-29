#ifndef WEMPCG_pwl
#define WEMPCG_pwl
/**************
 *  WEMPCG_pwl.h  *
 **************/


/*=============================================================*
 *  WEMPCG_pwl(A,b,x,epsi,W,F,p,M)                                 *
 *	                                                       *
 *  Verfahren der konjugierten Gradienten zur Loesung des      *
 *  linearen Gleichungssystems                                 *
 *	                                                       *
 *	      A1*x = -(A2^(-1)/(epsilon-1)+I)*b.               *
 *	                                                       *
 *  Vorkonditionierung per Wavelet-Preconditioner.             *
 *	                                                       *
 *  Parameter :                                                *
 *		A    : Matrix im sparse2-Format                *
 *		G    : Massenmatrix                            *
 *		b    : rechte Seite                            *
 *		x    : Startwert und Endwert                   *
 *		epsi : Genauigkeit                             *
 *		W    : Liste der Wavelets                      *
 *		F    : Elementliste der Einskalenbasis         *
 *		p    : Zahl der Paramtergebiete                *
 *		M    : Zahl der Level                          *
 *=============================================================*/


unsigned int WEMPCG_pwl(sparse2 *A, 
                        double *b, 
                        double *x, 
                        double epsi, 
                        wavelet_pwl *W, 
                        unsigned int **F, 
                        unsigned int p, 
                        unsigned int M);
#endif
