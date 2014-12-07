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
