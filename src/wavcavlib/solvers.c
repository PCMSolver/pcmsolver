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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"
 

void kahf_matr_lujs(double **MAT,int n,double *x,double *y)
{int i,k;
double temp;
for(i=0;i<n;i++)
	{temp=0.0;
	for(k=0;k<n;k++)
		temp=temp+MAT[i][k]*x[k];
	y[i]=temp;
	}
}



void nepv_test_howt(double **MAT,int n,double *x,double *rhs)
{int i;
double *y,err;
y=(double *)malloc(n*sizeof(double));
kahf_matr_lujs(MAT,n,x,y);

	
err=0.0;
for(i=0;i<n;i++)
	err=err+(rhs[i]-y[i]);
err=err/(double)n;
free(y);
fprintf(tmpout,"error=%f\n",err);
}


void vold_swap_puqb(int k,int j,double **A,double *b,int n)
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


int javs_find_gokm_pedt(double **MAT,double *b,double *x,int n)
{int i,j,k,h,m,p,n1,ind,res;
double l,pre,eps=1.0e-7,**A;
A=(double **)malloc(n*sizeof(double*));
for(i=0;i<n;i++)
	A[i]=(double *)malloc(n*sizeof(double));
for(i=0;i<n;i++)
	for(j=0;j<n;j++)
		A[i][j]=MAT[i][j];
for(m=0;m<n;m++)
   x[m]=0.0;

n1=n-1;
for(k=0;k<n1;k++)
	{h=k+1;
	ind=1;
	if(fabs(A[k][k])<eps)
		{p=k+1;
		while(p<n)
			{if(fabs(A[p][k])>=eps)
				{vold_swap_puqb(p,k,A,b,n);
				ind=2;
				}
			p++;
			}
		}    
	else
		ind=2;
	if(ind==2)
		{for(i=h;i<=n1;i++)
			{l=A[i][k]/A[k][k];
			for(j=k;j<=n1;j++)
				{A[i][j]=A[i][j]-(l*A[k][j]);
				}
			b[i]=b[i]-(l*b[k]);
			}
		}
	}

res=1;
for(i=0;i<n;i++)
if(fabs(A[i][i])<eps)
	{res=0;
	break;
	}

for(j=0;j<=n1;j++)
	{k=n1-j;
	x[k]=b[k];
	h=k+1;
	for(i=h;i<=n1;i++)
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
return res;
}


double fewl_diff_gurz(double *x,double *y,int n)
{int i;
double res;
res=0.0;
for(i=0;i<n;i++)
	res=res+fabs(x[i]-y[i]);
res=res/(double)n;
return res;
}



void dusr_find_vifq_dovz(double **A,int n,double *rhs,double *sol,int max,double acc)
{int i,j,k;
double *cur,*newit,temp,err;
cur=(double *)malloc(n*sizeof(double));
newit=(double *)malloc(n*sizeof(double));

for(i=0;i<n;i++)
	if(fabs(A[i][i])>0.0001)	cur[i]=1.0/A[i][i];
	else	cur[i]=1.0;

for(k=0;k<max;k++)
	{for(i=0;i<n;i++)
		{temp=0.0;
		for(j=0;j<n;j++)if(j!=i)
			temp=temp-A[i][j]*cur[j];
		temp=temp+rhs[i];
		newit[i]=temp/A[i][i];
		}
	
	err=fewl_diff_gurz(cur,newit,n);
	if(err<acc)
		{fprintf(tmpout,"Accuracy attained\n");
		break;
		}
	
	for(i=0;i<n;i++)
		cur[i]=newit[i];
	}
for(i=0;i<n;i++)
	sol[i]=newit[i];
free(cur);
free(newit);
}



void CG(double **A,double *rhs,int n,double *x,
double accuracy,int max_iter)
{int j,i,k;
double *p,*r,*temp,alpha,beta,aux2,aux1,error,aux;
p=(double*)malloc(n*sizeof(double));
r=(double*)malloc(n*sizeof(double));
temp=(double*)malloc(n*sizeof(double));



for(i=0;i<n;i++)
	x[i]=1.0;
for(i=0;i<n;i++)
	{aux=0.0;
	for(j=0;j<n;j++)
		aux=aux+A[i][j]*x[j];
	r[i]=rhs[i]-aux;
	}
for(i=0;i<n;i++)
	p[i]=r[i];



for(j=0;j<max_iter;j++)
	{aux=0.0;
	for(i=0;i<n;i++)
		aux=aux+(r[i]*r[i]);
	for(i=0;i<n;i++)
		{aux1=0.0;
		for(k=0;k<n;k++)
			aux1=aux1+(A[i][k]*p[k]);
		temp[i]=aux1;
		}
	aux1=0.0;
	for(i=0;i<n;i++)
		aux1=aux1+temp[i]*p[i];
	alpha=aux/aux1;
	for(i=0;i<n;i++)
		x[i]=x[i]+alpha*p[i];
	
	
	
	error=0.0;
	for(i=0;i<n;i++)
		{aux2=fabs(alpha*p[i]);
		error=error+(aux2/n);
		}
	if(error<accuracy)
		j=max_iter;
	for(i=0;i<n;i++)
		r[i]=r[i]-alpha*temp[i];
	aux1=0.0;
	for(i=0;i<n;i++)
		aux1=aux1+r[i]*r[i];
	beta=aux1/aux;
	for(i=0;i<n;i++)
		p[i]=r[i]+beta*p[i];
	}
free(p);
free(r);
free(temp);
}


