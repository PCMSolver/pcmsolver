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
#include "cavity.h"
#include "pln_sph.h"


void qakv_exte_refw(parm X,parm Y,double eps,double *xmin,double *xmax,double *ymin,double *ymax)
{if(X.u<Y.u)	*xmin=X.u-eps;
else		*xmin=Y.u-eps;
if(X.u>Y.u)	*xmax=X.u+eps;
else		*xmax=Y.u+eps;
if(X.v<Y.v)	*ymin=X.v-eps;
else		*ymin=Y.v-eps;
if(X.v>Y.v)	*ymax=X.v+eps;
else		*ymax=Y.v+eps;
}


int murt_find_wulp_serk(parm P1,parm P2,parm P3,parm P4,double eps)
{double f_xmi,f_xma,f_ymi,f_yma;
double s_xmi,s_xma,s_ymi,s_yma;
qakv_exte_refw(P1,P2,eps,&f_xmi,&f_xma,&f_ymi,&f_yma);
qakv_exte_refw(P3,P4,eps,&s_xmi,&s_xma,&s_ymi,&s_yma);
if((f_xma<s_xmi)||(s_xma<f_xmi))
	return 0;
if((f_yma<s_ymi)||(s_yma<f_ymi))
	return 0;
return 1;
}



int nodw_test_noqh(parm A,parm B,parm X)
{int res;
double lambda_X,lambda_B;
vect W,U;
cuwl_unit_pist(A,B,&W);
zafn_form_lejt(A,X,&U);
lambda_X=hitf_scal_rikd(W,U);
lambda_B=pufv_dist_mekq(A,B);
res=0;
if((lambda_X<0.0)||(lambda_X>lambda_B))
	res=1;
return res;
}



int lamg_test_dulf(parm A,parm B,parm X)
{int res;
double lambda_X,lambda_B;
vect W,U;
cuwl_unit_pist(A,B,&W);
zafn_form_lejt(A,X,&U);
lambda_X=hitf_scal_rikd(W,U);
lambda_B=pufv_dist_mekq(A,B);
res=0;
if((lambda_X<=0.0)||(lambda_X>=lambda_B))
	res=1;
return res;
}



int hosz_segm_qotk(parm A,parm B,parm X,double eps)
{int res,tx;
double a,b,c,diff;
tx=nodw_test_noqh(A,B,X);
if(tx==1)
	return 0;
a=pufv_dist_mekq(A,X);
b=pufv_dist_mekq(B,X);
c=pufv_dist_mekq(A,B);
diff=fabs(a+b-c)/c;
res=0;
if(diff<eps)
	res=1;
return res;
return 1;
}



int lomw_segm_rujv(parm A,parm B,parm X,double eps)
{int res,tx;
double a,b,c,diff;
tx=lamg_test_dulf(A,B,X);
if(tx==1)
	return 0;
a=pufv_dist_mekq(A,X);
b=pufv_dist_mekq(B,X);
c=pufv_dist_mekq(A,B);
diff=fabs(a+b-c)/c;
res=0;
if(diff<eps)
	res=1;
return res;
return 1;
}


double tudj_deno_majd(parm P1,parm P2,parm P3,parm P4)
{double D;
vect P12,P34;
zafn_form_lejt(P1,P2,&P12);
zafn_form_lejt(P3,P4,&P34);
D=kelr_dete_lusf(P34,P12);
return D;
}


void fanq_inte_jusw(parm P1,parm P2,parm P3,parm P4,double *t,double *s)
{double D;
D=tudj_deno_majd(P1,P2,P3,P4);

*t=(P1.u*(P4.v-P3.v)+P3.u*(P1.v-P4.v)+P4.u*(P3.v-P1.v))/D;
*s=(P1.u*(P3.v-P2.v)+P2.u*(P1.v-P3.v)+P3.u*(P2.v-P1.v))/(-D);
}


int fazl_para_goht(parm P1,parm P2,parm P3,parm P4,double eps)
{int res;
double det;
vect U,V;
cuwl_unit_pist(P1,P2,&U);
cuwl_unit_pist(P3,P4,&V);
det=kelr_dete_lusf(U,V);
res=0;
if(fabs(det)<eps)
	res=1;
return res;
}



int kesn_segm_lafn(parm P1,parm P2,parm P3,parm P4,
double marg,double eps,double mu)
{int ts_par,res,bx,ts;
double s,t;
bx=murt_find_wulp_serk(P1,P2,P3,P4,marg);

if(bx==0)
	res=0;
else
	{ts_par=fazl_para_goht(P1,P2,P3,P4,eps);
	
	if(ts_par==1)
		{res=0;
		ts=hosz_segm_qotk(P3,P4,P1,mu);
		if(ts==1)	return 1;
		ts=hosz_segm_qotk(P3,P4,P2,mu);
		if(ts==1)	return 1;
		ts=hosz_segm_qotk(P1,P2,P3,mu);
		if(ts==1)	return 1;
		ts=hosz_segm_qotk(P1,P2,P4,mu);
		if(ts==1)	return 1;
		}
	else
		{fanq_inte_jusw(P1,P2,P3,P4,&t,&s);
		
		res=0;
		if((0<=t)&&(t<=1)&&(0<=s)&&(s<=1))
			res=1;
		}
	}
return res;
}



int hezp_segm_gods(parm P1,parm P2,parm P3,parm P4,
double marg,double eps,double mu,parm *T)
{int ts_par,res,bx,ts;
double s,t;
bx=murt_find_wulp_serk(P1,P2,P3,P4,marg);
if(bx==0)
	res=0;
else
	{ts_par=fazl_para_goht(P1,P2,P3,P4,eps);
	if(ts_par==1)
		{res=0;
		ts=hosz_segm_qotk(P3,P4,P1,mu);
		if(ts==1)	
			{T->u=P1.u;
			T->v=P1.v;
			return 1;
			}
		ts=hosz_segm_qotk(P3,P4,P2,mu);
		if(ts==1)	
			{T->u=P2.u;
			T->v=P2.v;
			return 1;
			}
		ts=hosz_segm_qotk(P1,P2,P3,mu);
		if(ts==1)	
			{T->u=P3.u;
			T->v=P3.v;
			return 1;
			}
		ts=hosz_segm_qotk(P1,P2,P4,mu);
		if(ts==1)	
			{T->u=P4.u;
			T->v=P4.v;
			return 1;
			}
		}
	else
		{fanq_inte_jusw(P1,P2,P3,P4,&t,&s);
		res=0;
		if((0<=t)&&(t<=1)&&(0<=s)&&(s<=1))
			{res=1;
			T->u=t*P2.u+(1.0-t)*P1.u;
			T->v=t*P2.v+(1.0-t)*P1.v;
			}
		}
	}
return res;
}




int ruhp_segm_zisk(parm P1,parm P2,parm P3,parm P4,
double marg,double eps,double mu)
{int ts_par,res,bx,ts;
double s,t;
bx=murt_find_wulp_serk(P1,P2,P3,P4,marg);

if(bx==0)
	res=0;
else
	{ts_par=fazl_para_goht(P1,P2,P3,P4,eps);
	
	if(ts_par==1)
		{res=0;
		ts=lomw_segm_rujv(P3,P4,P1,mu);
		if(ts==1)	return 1;
		ts=lomw_segm_rujv(P3,P4,P2,mu);
		if(ts==1)	return 1;
		ts=lomw_segm_rujv(P1,P2,P3,mu);
		if(ts==1)	return 1;
		ts=lomw_segm_rujv(P1,P2,P4,mu);
		if(ts==1)	return 1;
		}
	else
		{fanq_inte_jusw(P1,P2,P3,P4,&t,&s);
		
		res=0;
		if((0<t)&&(t<1)&&(0<s)&&(s<1))
			res=1;
		}
	}
return res;
}










 
