/*****************
 *  WEMPGMRES.c  *
 *****************/


/*==============================================================*
 *  WEMPGMRES(A,b,x,epsi,p,M)                                   *
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
 *		p    : Anzahl der Paramtergebiete               *
 *		M    : Zahl der Level                           *
 *==============================================================*/


#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sparse2.h"
#include "precond.h"
#include "WEMPGMRES.h"


const unsigned int 	maxiter = 100;		/* Restart - Parameter */


/*=================*
 *  Hauptprogramm  *
 *=================*/
 
unsigned int WEMPGMRES1(A,b,x,epsi,p,M)
/* loest A2'*x = b */
sparse2 	*A;
double          *b, *x, epsi;
unsigned int	p, M;
{
signed int    	i, j;
unsigned int	k, l;
double  	rn, t;
double  	**H, **V, *D, *v, *f, *u, *c, *s, *y;

V = (double**) malloc((maxiter+1)*sizeof(double*));
H = (double**) malloc((maxiter+1)*sizeof(double*));

for (i=0; i<=maxiter; i++)  
{  V[i] = (double*) malloc(A->n*sizeof(double));
   H[i] = (double*) malloc(maxiter*sizeof(double));
   }
   
v = (double*) malloc(A->n*sizeof(double));
u = (double*) malloc((maxiter+1)*sizeof(double));
f = (double*) malloc((maxiter+1)*sizeof(double));
c = (double*) malloc(maxiter*sizeof(double));
s = (double*) malloc(maxiter*sizeof(double));
y = (double*) malloc(maxiter*sizeof(double));

/* Preconditioning: D = diag(A)^(-1) */
D = (double*) malloc(A->n*sizeof(double));
for (i=0; i<A->n; i++) 
{  get_sparse2(A,i,i,&rn,&t);
   D[i] = 1/t;
   }
rn = 1;
k = 0;

/* starte Iteration */
for (l=0; rn>epsi; l++)
{  
   /* Startresiduum: V(0,:) = D*(b - A*x), rn = sqrt(V(0,:)*V(0,:)') */
   rn = 0;
   memcpy(V[0],b,A->n*sizeof(double));
   for (i=0; i<A->n; i++)
   {  for (j=0; j<A->row_number[i]; j++)
      {  V[0][A->index[i][j]] -= A->value2[i][j] * x[i];
         }
      }
   for (i=0; i<A->n; i++)
   {  V[0][i] *= D[i];
      rn += V[0][i]*V[0][i];
      }
   rn = sqrt(rn);

   /* V(0,:) /= rn */
   t = 1/rn;
   for (i=0; i<A->n; i++) V[0][i] *= t;
   u[0] = rn;

   /* innere Iteration */
   for (k=0; (k<maxiter) && (rn>epsi); k++)
   {  
      /* Matrix-Vektor-Produkt v = D*(A*V(k,:)) */
      memset(v,0,A->n*sizeof(double));
      for (i=0; i<A->n; i++)
      {  for (j=0; j<A->row_number[i]; j++)
         {  v[A->index[i][j]] += A->value2[i][j] * V[k][i];
	    }
	 }
      for (i=0; i<A->n; i++) v[i] *= D[i];	 

      /* Berechnung der Skalps f = V(0:k,:)*v */
      memset(f,0,(k+1)*sizeof(double));
      for (i=0; i<=k; i++)
      {  for (j=0; j<A->n; j++) f[i] += V[i][j] * v[j];
	 }

      /* (k+1)-ter Basisvektor v = v - V(0:k,:)*f und rn = sqrt(v'*v) */
      rn = 0;
      for (i=0; i<A->n; i++)
      {  for (j=0; j<=k; j++) v[i] -= V[j][i] * f[j];
	 rn += v[i]*v[i];
         }
      rn = sqrt(rn);

      /* V(k+1,:) = v' / rn */
      t = 1/rn;
      for (i=0; i<A->n; i++) V[k+1][i] = v[i] * t;

      /* Berechnung der Givensrotationen */
      t = f[0];
      for (i=0; i<k; i++)
      {  f[i] = c[i] * t + s[i] * f[i+1];
	 t = c[i] * f[i+1] - s[i] * t;
	 }
      f[k] = sqrt(t*t + rn*rn);
      c[k] = t  / f[k];
      s[k] = rn / f[k];

      /* H(0:k,k) = f */
      for (i=0; i<=k; i++) H[i][k] = f[i];

      u[k+1] = - s[k] * u[k];
      u[k]   = c[k] * u[k];
      rn = fabs(u[k+1]);
      }   /* Ende der inneren Iteration */

   /* Loesen des Dreiecksystems y = H(0:k-1,0:k-1) \ u(0:k-1) */
   for (i=k-1; i>=0; i--)
   {  t = 0;
      for (j=k-1; j>i; j--) t += H[i][j] * y[j];
      y[i] = (u[i] - t) / H[i][i];
      }

   /* neue Startnaeherung x = x + V(0:k-1,:)'*y */
   for (i=0; i<A->n; i++)
   {  for (j=0; j<k; j++) x[i] += V[j][i] * y[j];
      }
   }

/* Speicherplatz wieder freigeben */
for (i=0; i<=maxiter; i++)
{  free(V[i]);
   free(H[i]);
   }

free(V);
free(H);
free(D);
free(v);
free(f);
free(u); 
free(c);
free(s);
free(y);

return((l-1)*maxiter+k);
}


unsigned int WEMPGMRES2(A,b,x,epsi,p,M)
/* loest A2*x = b */
sparse2 	*A;
double          *b, *x, epsi;
unsigned int	p, M;
{
signed int    	i, j;
unsigned int	k, l;
double  	rn, t;
double  	**H, **V, *D, *v, *f, *u, *c, *s, *y;

V = (double**) malloc((maxiter+1)*sizeof(double*));
H = (double**) malloc((maxiter+1)*sizeof(double*));

for (i=0; i<=maxiter; i++)  
{  V[i] = (double*) malloc(A->n*sizeof(double));
   H[i] = (double*) malloc(maxiter*sizeof(double));
   }
   
v = (double*) malloc(A->n*sizeof(double));
u = (double*) malloc((maxiter+1)*sizeof(double));
f = (double*) malloc((maxiter+1)*sizeof(double));
c = (double*) malloc(maxiter*sizeof(double));
s = (double*) malloc(maxiter*sizeof(double));
y = (double*) malloc(maxiter*sizeof(double));

/* Preconditioning: D = diag(A)^(-1) */
D = (double*) malloc(A->n*sizeof(double));
for (i=0; i<A->n; i++) 
{  get_sparse2(A,i,i,&rn,&t);
   D[i] = 1/t;
   }
rn = 1;
k = 0;

/* starte Iteration */
for (l=0; rn>epsi; l++)
{  
   /* Startresiduum: V(0,:) = D*(b - A*x), rn = sqrt(V(0,:)*V(0,:)') */
   rn = 0;
   memcpy(V[0],b,A->n*sizeof(double));
   for (i=0; i<A->n; i++)
   {  for (j=0; j<A->row_number[i]; j++)
      {  V[0][i] -= A->value2[i][j] * x[A->index[i][j]];
         }
      V[0][i] *= D[i];
      rn += V[0][i]*V[0][i];
      }
   rn = sqrt(rn);

   /* V(0,:) /= rn */
   t = 1/rn;
   for (i=0; i<A->n; i++) V[0][i] *= t;
   u[0] = rn;

   /* innere Iteration */
   for (k=0; (k<maxiter) && (rn>epsi); k++)
   {  
      /* Matrix-Vektor-Produkt v = D*(A*V(k,:)) */
      memset(v,0,A->n*sizeof(double));
      for (i=0; i<A->n; i++)
      {  for (j=0; j<A->row_number[i]; j++)
         {  v[i] += A->value2[i][j] * V[k][A->index[i][j]];
	    }
	 v[i] *= D[i];
	 }

      /* Berechnung der Skalps f = V(0:k,:)*v */
      memset(f,0,(k+1)*sizeof(double));
      for (i=0; i<=k; i++)
      {  for (j=0; j<A->n; j++) f[i] += V[i][j] * v[j];
	 }

      /* (k+1)-ter Basisvektor v = v - V(0:k,:)*f und rn = sqrt(v'*D*v) */
      rn = 0;
      for (i=0; i<A->n; i++)
      {  for (j=0; j<=k; j++) v[i] -= V[j][i] * f[j];
	 rn += v[i]*v[i];
         }
      rn = sqrt(rn);

      /* V(k+1,:) = v' / rn */
      t = 1/rn;
      for (i=0; i<A->n; i++) V[k+1][i] = v[i] * t;

      /* Berechnung der Givensrotationen */
      t = f[0];
      for (i=0; i<k; i++)
      {  f[i] = c[i] * t + s[i] * f[i+1];
	 t = c[i] * f[i+1] - s[i] * t;
	 }
      f[k] = sqrt(t*t + rn*rn);
      c[k] = t  / f[k];
      s[k] = rn / f[k];

      /* H(0:k,k) = f */
      for (i=0; i<=k; i++) H[i][k] = f[i];

      u[k+1] = - s[k] * u[k];
      u[k]   = c[k] * u[k];
      rn = fabs(u[k+1]);
      }   /* Ende der inneren Iteration */

   /* Loesen des Dreiecksystems y = H(0:k-1,0:k-1) \ u(0:k-1) */
   for (i=k-1; i>=0; i--)
   {  t = 0;
      for (j=k-1; j>i; j--) t += H[i][j] * y[j];
      y[i] = (u[i] - t) / H[i][i];
      }

   /* neue Startnaeherung x = x + V(0:k-1,:)'*y */
   for (i=0; i<A->n; i++)
   {  for (j=0; j<k; j++) x[i] += V[j][i] * y[j];
      }
   }

/* Speicherplatz wieder freigeben */
for (i=0; i<=maxiter; i++)
{  free(V[i]);
   free(H[i]);
   }

free(V);
free(H);
free(D);
free(v);
free(f);
free(u); 
free(c);
free(s);
free(y);

return((l-1)*maxiter+k);
}


/*==============================================================*
 *  WEMPGMRES(A,B,rhs,x,epsi,p,M)                               *
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
 *		p    : Anzahl der Paramtergebiete               *
 *		M    : Zahl der Level                           *
 *==============================================================*/

unsigned int WEMPGMRES3(A,B,rhs,x,epsi,p,M)
/* loest (B1*G^(-1)*A2'-B2*G^(-1)*A1)*x = rhs */
sparse2 	*A, *B;
double          *rhs, *x, epsi;
unsigned int	p, M;
{
signed int    	i, j;
unsigned int	k, l;
double  	rn, t;
double  	**H, **V, *v, *f, *u, *c, *s, *y;

/* Speicherplatz allokieren */
V = (double**) malloc((maxiter+1)*sizeof(double*));
H = (double**) malloc((maxiter+1)*sizeof(double*));

for (i=0; i<=maxiter; i++)  
{  V[i] = (double*) malloc(A->n*sizeof(double));
   H[i] = (double*) malloc(maxiter*sizeof(double));
   }
   
v = (double*) malloc(A->n*sizeof(double));
u = (double*) malloc((maxiter+1)*sizeof(double));
f = (double*) malloc((maxiter+1)*sizeof(double));
c = (double*) malloc(maxiter*sizeof(double));
s = (double*) malloc(maxiter*sizeof(double));
y = (double*) malloc(maxiter*sizeof(double));

rn = 1;
k = 0;

/* starte Iteration */
for (l=0; rn>epsi; l++)
{  
   /* Startresiduum: V(0,:) = rhs - (B1*G^(-1)*A2'-B2*G^(-1)*A1)*x */
   memset(v,0,A->n*sizeof(double));
   memcpy(V[0],rhs,A->n*sizeof(double));
   for (i=0; i<A->n; i++)
   {  for (j=0; j<A->row_number[i]; j++)
      {  v[A->index[i][j]] += A->value2[i][j] * x[i];
         }
      }
   inv_G_times_x(v,p,M);
   for (i=0; i<B->n; i++)
   {  for (j=0; j<B->row_number[i]; j++)
      {  V[0][i] -= B->value1[i][j] * v[B->index[i][j]];
         }
      }
   memset(v,0,A->n*sizeof(double));
   for (i=0; i<A->n; i++)
   {  for (j=0; j<A->row_number[i]; j++)
      {  v[i] += A->value1[i][j] * x[A->index[i][j]];
         }
      }
   inv_G_times_x(v,p,M);
   for (i=0; i<B->n; i++)
   {  for (j=0; j<B->row_number[i]; j++)
      {  V[0][i] += B->value2[i][j] * v[B->index[i][j]];
         }
      }	 
   	 
   /* v = precond(V(0,:)) und rn = sqrt(v'*v) */ 	 
   rn = 0;
   precond(v,V[0],p,M);
   for (i=0; i<A->n; i++) rn += v[i]*v[i];
   rn = sqrt(rn);

   /* v = V(0,:)/rn */
   t = 1/rn;
   for (i=0; i<A->n; i++) V[0][i] = t * v[i];
   u[0] = rn;

   /* innere Iteration */
   for (k=0; (k<maxiter) && (rn>epsi); k++)
   {  
      /* Matrix-Vektor-Produkt V[k+1][i] = (B1*G^(-1)*A2'-B2*G^(-1)*A1)*V(k,:) */
      memset(v,0,A->n*sizeof(double));
      memset(V[k+1],0,A->n*sizeof(double));
      for (i=0; i<A->n; i++)
      {  for (j=0; j<A->row_number[i]; j++)
         {  v[A->index[i][j]] += A->value2[i][j] * V[k][i];
            }
         }
      inv_G_times_x(v,p,M);
      for (i=0; i<B->n; i++)
      {  for (j=0; j<B->row_number[i]; j++)
         {  V[k+1][i] += B->value1[i][j] * v[B->index[i][j]];
            }
         }
      memset(v,0,A->n*sizeof(double));
      for (i=0; i<A->n; i++)
      {  for (j=0; j<A->row_number[i]; j++)
         {  v[i] += A->value1[i][j] * V[k][A->index[i][j]];
            }
         }
      inv_G_times_x(v,p,M);
      for (i=0; i<B->n; i++)
      {  for (j=0; j<B->row_number[i]; j++)
         {  V[k+1][i] -= B->value2[i][j] * v[B->index[i][j]];
            }
         }

      /* v = precond(V(k+1,:) */
      precond(v,V[k+1],p,M);
      
      /* Berechnung der Skalps f = V(0:k,:)*v */
      memset(f,0,(k+1)*sizeof(double));
      for (i=0; i<=k; i++)
      {  for (j=0; j<A->n; j++) f[i] += V[i][j] * v[j];
	 }

      /* (k+1)-ter Basisvektor v = v - V(0:k,:)*f und rn = sqrt(v'*v) */
      rn = 0;
      for (i=0; i<A->n; i++)
      {  for (j=0; j<=k; j++) v[i] -= V[j][i] * f[j];
	 rn += v[i]*v[i];
         }
      rn = sqrt(rn);

      /* V(k+1,:) = v' / rn */
      t = 1/rn;
      for (i=0; i<A->n; i++) V[k+1][i] = v[i] * t;

      /* Berechnung der Givensrotationen */
      t = f[0];
      for (i=0; i<k; i++)
      {  f[i] = c[i] * t + s[i] * f[i+1];
	 t = c[i] * f[i+1] - s[i] * t;
	 }
      f[k] = sqrt(t*t + rn*rn);
      c[k] = t  / f[k];
      s[k] = rn / f[k];

      /* H(0:k,k) = f */
      for (i=0; i<=k; i++) H[i][k] = f[i];

      u[k+1] = - s[k] * u[k];
      u[k]   = c[k] * u[k];
      rn = fabs(u[k+1]);
      }   /* Ende der inneren Iteration */

   /* Loesen des Dreiecksystems y = H(0:k-1,0:k-1) \ u(0:k-1) */
   for (i=k-1; i>=0; i--)
   {  t = 0;
      for (j=k-1; j>i; j--) t += H[i][j] * y[j];
      y[i] = (u[i] - t) / H[i][i];
      }

   /* neue Startnaeherung x = x + V(0:k-1,:)'*y */
   for (i=0; i<A->n; i++)
   {  for (j=0; j<k; j++) x[i] += V[j][i] * y[j];
      }
   }

/* Speicherplatz wieder freigeben */
for (i=0; i<=maxiter; i++)
{  free(V[i]);
   free(H[i]);
   }
free(V);
free(H);
free(v);
free(f);
free(u); 
free(c);
free(s);
free(y);
return((l-1)*maxiter+k);
}
