/*
 * Purpose: Patch representation of molecular cavities from
 *			atomic coordinates and radii.
 * Version: September 21, 2009.
 * Author : Maharavo Randrianarivony.
 * Affiliation: Institut für Numerische Simulation.
 *              Rheinische Friedrich-Wilhelm Universität Bonn.
 *              Wegelerstraße 6, Bonn 53115.
 *              Germany.
 */

#include <math.h>
#include <stdlib.h>
#include "cavity.h"
#include "splinemol.h"


void behn_swap_fujk(int k,int j,double **A,double *b,int n)
{int i;
double temp;
for(i=0;i<n;i++)
	{temp=A[k][i];
	A[k][i]=A[j][i];
	A[j][i]=temp;
	}
temp=b[k];
b[k]=b[j];
b[j]=temp;
}

 
void poml_trid_fert(double **MAT,double *rhs,double *x,int n)
{int i,j,k,h,m,p,n1,ind,mx;
double l,pre,eps=1.0e-7,**A,*b;

A=(double **)malloc(n*sizeof(double*));
for(i=0;i<n;i++)
	A[i]=(double *)malloc(n*sizeof(double));
b=(double *)malloc(n*sizeof(double));
for(i=0;i<n;i++)
	{for(j=0;j<n;j++)
		A[i][j]=MAT[i][j];
	b[i]=rhs[i];
	}
for(m=0;m<n;m++)
   x[m]=0.0;

n1=n-1;
for(k=0;k<n1;k++)
	{ind=1;
	if(fabs(A[k][k])<eps)
		{p=k+1;
		while(p<n)
			{if(fabs(A[p][k])>=eps)
				{behn_swap_fujk(p,k,A,b,n);
				ind=2;
				}
			p++;
			}
		}    
	else
		ind=2;
	if(ind==2)
		{l=A[k+1][k]/A[k][k];
		mx=k+2;
		if(mx>n1)	mx=n1;
		for(j=k;j<=mx;j++)
			A[k+1][j]=A[k+1][j]-(l*A[k][j]);
		b[k+1]=b[k+1]-(l*b[k]);
		}
	}

for(j=0;j<=n1;j++)
	{k=n1-j;
	x[k]=b[k];
	h=k+1;
	mx=h+2;
	if(mx>n1)	mx=n1;
	for(i=h;i<=mx;i++)
		{pre=A[k][i]*x[i];
		x[k]=x[k]-pre;
		}
	if(fabs(A[k][k])<eps)
		x[k]=0.1;
	else
		x[k]=x[k]/A[k][k];
	}

for(i=0;i<n;i++)
	free(A[i]);
free(A);
free(b);
}
