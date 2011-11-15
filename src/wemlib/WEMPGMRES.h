#ifndef WEMPGMRES
#define WEMPGMRES
/*****************
 *  WEMPGMRES.h  *
 *****************/


/*==============================================================*
 *  WEMPGMRES(A,b,x,epsi,p,M)                                   *
 *	                                                        *
 *  GMRES-Verfahren zur Loesung des linearen Gleichungssystems  *
 *	                                                        *
 *		A2'*x = b  bzw. A2*x = b.                       *
 *	                                                        *
 *  Vorkonditionierung per Diagonalskalierung.                  *
 *	                                                        *
 *  Parameter :                                                 *
 *		A    : Matrix im sparse2-Format                 *
 *		b    : rechte Seite                             *
 *		x    : Startwert und Endwert                    *
 *		epsi : Genauigkeit                              *
 *		p    : Anzahl der Paramtergebiete               *
 *		M    : Zahl der Level                           *
*===============================================================*/


unsigned int WEMPGMRES1(sparse2 *A, double *b, double *x, double epsi, unsigned int p, unsigned int M);


unsigned int WEMPGMRES2(sparse2 *A, double *b, double *x, double epsi, unsigned int p, unsigned int M);


/*==============================================================*
 *  WEMPGMRES(A,B,rhs,x,epsi,p,M)                               *
 *	                                                        *
 *  GMRES-Verfahren zur Loesung des linearen Gleichungssystems  *
 *	                                                        *
 *	    (B1*G^(-1)*A2^T-B2*G^(-1)*A1)*x = rhs.              *
 *	                                                        *
 *  Vorkonditionierung per Wavelet-Preconditioner.              *
 *	                                                        *
 *  Parameter :                                                 *
 *		A, B : Matrizen im sparse2-Format               *
 *		rhs  : rechte Seite                             *
 *		x    : Startwert und Endwert                    *
 *		epsi : Genauigkeit                              *
 *		p    : Anzahl der Paramtergebiete               *
 *		M    : Zahl der Level                           *
 *==============================================================*/

unsigned int WEMPGMRES3(sparse2 *A, sparse2 *B, double *rhs, double *x, double epsi, unsigned int p, unsigned int M);
#endif
