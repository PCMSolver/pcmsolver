#include <iostream>
#include <ostream>
#include <cmath>
#include <cstdlib>

extern "C" {
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_legendre.h>
}

#include <Eigen/Dense>

using namespace Eigen;

using namespace std;

#include "Getkw.h"
#include "GreensFunction.h"
#include "DielectricSphere.h"

int    DielectricSphere::nSteps = 1000;
int    DielectricSphere::maxLGreen = 30;    //max L in getting the final summation of Green's function
int    DielectricSphere::maxLEpsilon = 60;  //max L in getting the constant of C(r1,r2)
double DielectricSphere::hStep = 5.0e-4;    //step size for integration
double DielectricSphere::rBegin = 0.5;   //start poinr for solving the differential equation in the radial dimension
double DielectricSphere::rEnd = 103.0;    //final point for solving the differential equation in the radial dimension
double DielectricSphere::rMin = 1.0;     //min scanning point in radial part in solving differential equation
double DielectricSphere::rMax = 100.0;    //max scanning point in radial part in solving differential equation

typedef struct{
  DielectricSphere *cont;
  int L;
} GSL_param_helper;

int evalFunc(double t, const double ya[], double f[], void *params){  
	GSL_param_helper *help = (GSL_param_helper *)params;
	DielectricSphere *cont = help->cont;
	int L = help->L;
	double er, der_er;
	cont->profile(&er, &der_er, t);
	
	f[0] = ya[1];                                               // y' = z
	f[1] = L*(L+1.0)/(t*t)-ya[1]*ya[1]-(2.0/t+der_er/er)*ya[1]; // y'' = l(l+1)/t^2  -(2/t+der_er/er+y')y'--derivation by Cao
	
	return GSL_SUCCESS;
}

DielectricSphere::DielectricSphere(double epsIn, double epsOut, 
								   Vector3d posSph, double radSph, 
								   double widthInt)
{
	epsInside = epsIn;
	epsOutside = epsOut;
	spherePosition = posSph;
	sphereRadius = radSph;
	interfaceWidth = widthInt;
	initDielectricSphere();
};

DielectricSphere::DielectricSphere(Section green)
{
	epsInside = green.getDbl("EpsIn");
	epsOutside = green.getDbl("EpsOut");
	sphereRadius = green.getDbl("SphereRadius");
	interfaceWidth = green.getDbl("InterfaceWidth");
	const vector<double> &pos_ = green.getDblVec("SpherePosition");
	spherePosition << pos_[0], pos_[1], pos_[2];
	initDielectricSphere();
}

void DielectricSphere::initDielectricSphere() {
	double deltaR = rEnd - rBegin;
	nSteps = deltaR / hStep + 1;
	grid = VectorXd::Zero(nSteps);

	for(int i = 0; i <= maxLGreen; i++) {
		FuncGrid f1 = FuncGrid::Zero(2, nSteps);
		FuncGrid f2 = FuncGrid::Zero(2, nSteps);
		radialG1.push_back(f1);
		radialG2.push_back(f2);
	}
	radialC1 = FuncGrid::Zero(2, nSteps);
	radialC2 = FuncGrid::Zero(2, nSteps);

	for (int i = 0; i < nSteps; i++) {
		grid(i) = rBegin + i * hStep;
	}
	rEnd = grid(nSteps - 1);
	independent_solutions();
	uniformFlag = false;
}

double DielectricSphere::evalf(Vector3d &p1, Vector3d &p2) {
	return getSingleLayer(p1, p2);
}

double DielectricSphere::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2) {
	exit(-1);
	return 0.0;
	//	return getDoubleLayer(p1, p2, direction);
}


double DielectricSphere::getSingleLayer(Vector3d &p1, Vector3d &p2) {
	double r1 = p1.norm();
	double r2 = p2.norm();
	double r12= (p1-p2).norm();
	double cos_theta = p1.dot(p2)/(r1*r2);
	double gfsph = converged_gf(r1, r2, r12, cos_theta); 
	return gfsph;
}

void DielectricSphere::independent_solutions(){

	cout << "Computing coeffiencents..." << endl;
	getU1(rEnd, rBegin, radialC1, maxLEpsilon);
	getU2(rBegin, rEnd, radialC2, maxLEpsilon);

	cout << "Computing first radial solution..." << endl;
	for (int l = 0; l <= maxLGreen; l++){
		getU1(rEnd, rBegin, radialG1[l], l);
	}

	cout << "Computing second radial solution..." << endl;
	for (int l = 0; l <= maxLGreen; l++){
		getU2(rBegin, rEnd, radialG2[l], l);
	}
}

void DielectricSphere::getU1(double r, double r0, FuncGrid & points, int L){

	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step * s  = gsl_odeiv_step_alloc (T, 2);
	GSL_param_helper help;
	help.cont = this;
	help.L = L;
	gsl_odeiv_system sys = {evalFunc, NULL, 2, (void *) &help};  //the fourth parameter to be this object
	double t = r0, t1 = r;
	double ya[2];

	ya[0] = L*log(r0);   // y
	ya[1] = L/r0; // y'
	double y_err[2];
	double dydt_in[2], dydt_out[2];
	//	printf("start %15.10e %15.10e %15.10e\n", t, ya[0], ya[1]);
	GSL_ODEIV_FN_EVAL(&sys, t, ya, dydt_in);
	for(int i = 0; i < nSteps; i++) { 
		points(0, i) = ya[0];
		points(1, i) = ya[1];
		t = grid(i);
		int status = gsl_odeiv_step_apply (s, t, hStep, ya, y_err, dydt_in,
										   dydt_out,&sys);
		dydt_in[0] = dydt_out[0];
		dydt_in[1] = dydt_out[1];
		//		printf("steps %15.10e %15.10e %15.10e\n", grid[i], points(0,i), points(1,i));
		if (status != GSL_SUCCESS) {
			cout << "getU1 error" << endl;
			exit(-1);
		}
	}
	gsl_odeiv_step_free(s);
}

void DielectricSphere::getU2(double r, double r0, FuncGrid & points, int L){

	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step * s  = gsl_odeiv_step_alloc (T, 2);
	GSL_param_helper help;
	help.cont = this;
	help.L = L;
	gsl_odeiv_system sys = {evalFunc, NULL, 2, &help};  //the fourth parameter to be this object
	double t = r0, t1 = r;
	double ya[2];

	ya[0] = -(L+1)*log(r0);   // y
	ya[1] = -(L+1)/r0; // y'
	double y_err[2];
	double dydt_in[2], dydt_out[2];
	GSL_ODEIV_FN_EVAL(&sys, t, ya, dydt_in);
	for(int i = nSteps - 1; i >= 0; i--) {
		points(0, i) = ya[0];
		points(1, i) = ya[1];
		t = grid(i);
		int status = gsl_odeiv_step_apply (s, t, -hStep, ya, y_err, dydt_in, 
										   dydt_out,&sys);
		dydt_in[0] = dydt_out[0];
		dydt_in[1] = dydt_out[1];
		if (status != GSL_SUCCESS) {
			cout << "getU2 error" << endl;
			exit(-1);
		}
	}
	gsl_odeiv_step_free(s);

}

double DielectricSphere::converged_gf(double r1, double r2, double r12, 
								  double cos_theta){
	double Cr12 = 0.0;
	double eps_r2, deps_r2;
	profile(&eps_r2, &deps_r2, r2);

	double gr12=greenfunc(r1, r2, r12, cos_theta, radialC1, radialC2, maxLEpsilon, eps_r2, 0, Cr12);

	if (r1<r2) {
		Cr12=pow(r1/r2,maxLEpsilon)/(gr12*r2);
	} else {
		Cr12=pow(r2/r1,maxLEpsilon)/(gr12*r1);
	}

	gr12=0.0;
	for (int l = 0; l < maxLGreen; l++) {
		gr12 += greenfunc(r1, r2, r12, cos_theta, radialG1[l], radialG2[l], l, eps_r2, 1, Cr12);
	}
	double gr = 1.0 / (Cr12 * r12) + gr12;
	return gr;
}

void DielectricSphere::profile(double *er, double *der_er, double r){
	double tanh_r = tanh((r-sphereRadius)/interfaceWidth);     //function tanh((r-cen_r0)/param_D)
	*er = 0.5*(epsInside + epsOutside) 
		+ 0.5*(epsOutside - epsInside) * tanh_r;     //epsilon(r)
	*der_er = 0.5*(epsInside - epsOutside) 
		* ( 1 - tanh_r * tanh_r) 
		/ interfaceWidth; //first derivative of epsilon(r)
}        

double DielectricSphere::greenfunc(double r1, double r2, double r12, 
								   double cos_theta, FuncGrid f1, FuncGrid f2, 
								   int l, double epi_r2, int param_Cr12, 
								   double Cr12){

	double g12 = 0.0;
	double u11, u12, u21, u22, d11, d12, d21, d22; 
	double plx = gsl_sf_legendre_Pl(l, cos_theta);

	int idx1 = int((r1-rBegin) / hStep) - 1;
	int idx2 = int((r2-rBegin) / hStep) - 1; 

	double delta1 = r1 - grid(idx1);
	double delta2 = r2 - grid(idx2);

	/*
	cout << "Gridpoints " << " " << grid(idx1)  << " " << grid(idx2) << " " << delta1 << " " << delta2 << endl;
	cout << "Indices " << idx1 << " " << idx2 << endl;
	cout << "Functions " << grid(idx1+1) << " " << f1(0,idx1+1) << " " << f1(1,idx1+1) << endl; // i = 0
	cout << "Functions " << grid(idx1+1) << " " << f2(0,idx1+1) << " " << f2(1,idx1+1) << endl;
	cout << "Functions " << grid(idx2+1) << " " << f1(0,idx2+1) << " " << f1(1,idx2+1) << endl; //2
	cout << "Functions " << grid(idx2+1) << " " << f2(0,idx2+1) << " " << f2(1,idx2+1) << endl;
	*/
		
	u11 = f1(0, idx1) + (f1(0, idx1 + 1) - f1(0, idx1)) * delta1 / hStep;	    								  
	u12 = f1(0, idx2) + (f1(0, idx2 + 1) - f1(0, idx2)) * delta2 / hStep;	    								    
	u21 = f2(0, idx1) + (f2(0, idx1 + 1) - f2(0, idx1)) * delta1 / hStep;	    								    
	u22 = f2(0, idx2) + (f2(0, idx2 + 1) - f2(0, idx2)) * delta2 / hStep;	    								    

	d11 = f1(1, idx1) + (f1(1, idx1 + 1) - f1(1, idx1)) * delta1 / hStep;	    								    
	d12 = f1(1, idx2) + (f1(1, idx2 + 1) - f1(1, idx2)) * delta2 / hStep;
	d21 = f2(1, idx1) + (f2(1, idx1 + 1) - f2(1, idx1)) * delta1 / hStep;	    								    
	d22 = f2(1, idx2) + (f2(1, idx2 + 1) - f2(1, idx2)) * delta2 / hStep;

	if (r1<r2) {
		g12 = exp(u11-u12)*(2*l+1)/((d12-d22)*r2*r2*epi_r2);
		if (param_Cr12 == 1) g12 = (g12 - pow(r1/r2,l)/(r2*Cr12))*plx;
	} else {
		g12 = exp(u21-u22)*(2*l+1)/((d12-d22)*r2*r2*epi_r2);
		if (param_Cr12 == 1) g12 = (g12 - pow(r2/r1,l)/(r1*Cr12))*plx;
	}

	return g12;
}

double DielectricSphere::getDoubleLayer(Vector3d &p1, Vector3d &p2, Vector3d &n2){
	Vector3d grad = gradient(p1, p2);
	double eps2, deps2;
	profile(&eps2, &deps2, p2.norm());
	double dergfsph = eps2 * grad.dot(n2);
	return dergfsph;
}

void DielectricSphere::grad(Vector3d &g, Vector3d &p1, Vector3d &p2){
	double cos_theta = p1.dot(p2)/(p1.norm() * p2.norm());
	double * plx = new double[maxLEpsilon+1];     //for storing the lengedre function
	double * dplx = new double[maxLEpsilon+1];    //for storing the first derivative of lengedre function	
	gsl_sf_legendre_Pl_deriv_array(maxLEpsilon, cos_theta, plx, dplx);
	g = converged_deri_gf(p1, p2, plx, dplx);
	//	cout << "Gradient " << g.transpose() << endl;
	delete plx;
	delete dplx;
	return;
}

Vector3d DielectricSphere::converged_deri_gf(VectorXd p1, VectorXd p2, double *plx, 
											 double *dplx){
	Vector3d greenGradient(0, 0, 0);
	Vector3d Cr12Gradient(0, 0, 0);
	double r1 = p1.norm();
	double r2 = p2.norm();
	double eps2, deps2;
	profile(&eps2, &deps2, r2);

	double pl = plx[maxLEpsilon];
	double dpl = dplx[maxLEpsilon];
	
	Cr12Gradient.Zero();
	greenGradient = greenfunc_der(p1, p2, Cr12Gradient, radialC1, radialC2, 
								  pl, dpl, maxLEpsilon, 0);
	if (r1<r2) {
		double powl = pow(r1/r2, maxLEpsilon);
		Cr12Gradient(0)=powl*(pl*p2(0)*(-maxLEpsilon-1)+dpl*(p1(0)*r2*r2-p2(0)*(p1.dot(p2)))/(r1*r2))/(r2*r2*r2*greenGradient(0));
		Cr12Gradient(1)=powl*(pl*p2(1)*(-maxLEpsilon-1)+dpl*(p1(1)*r2*r2-p2(1)*(p1.dot(p2)))/(r1*r2))/(r2*r2*r2*greenGradient(1));
		Cr12Gradient(2)=powl*(pl*p2(2)*(-maxLEpsilon-1)+dpl*(p1(2)*r2*r2-p2(2)*(p1.dot(p2)))/(r1*r2))/(r2*r2*r2*greenGradient(2));
	} else {
		double powl = pow(r2/r1, maxLEpsilon);
		Cr12Gradient(0)=powl*(pl*p2(0)*maxLEpsilon+dpl*(p1(0)*r2*r2-p2(0)*(p1.dot(p2)))/(r1*r2))/(r1*r2*r2*greenGradient(0));
		Cr12Gradient(1)=powl*(pl*p2(1)*maxLEpsilon+dpl*(p1(1)*r2*r2-p2(1)*(p1.dot(p2)))/(r1*r2))/(r1*r2*r2*greenGradient(1));
		Cr12Gradient(2)=powl*(pl*p2(2)*maxLEpsilon+dpl*(p1(2)*r2*r2-p2(2)*(p1.dot(p2)))/(r1*r2))/(r1*r2*r2*greenGradient(2));
	}
	
	greenGradient.Zero();
	
	for (int l = 0; l < maxLGreen; l++) {
		Vector3d gl = greenfunc_der(p1, p2, Cr12Gradient, radialG1[l],
									radialG2[l], plx[l], dplx[l], l, 1);
		greenGradient += gl;
	}
	
	double r12 = (p1-p2).norm();
	Vector3d gr_d;
	gr_d.array() = (p1.array() - p2.array()) / 
		(Cr12Gradient.array() * r12 * r12 * r12) + greenGradient.array();
	gr_d *= -1;
	return gr_d;
}

Vector3d DielectricSphere::greenfunc_der(Vector3d p1, Vector3d p2, Vector3d Cr12,
										 FuncGrid f1, FuncGrid f2, double plx, 
										 double dplx, int l, int flagCr12)
{
	double eps2, deps2;
	double prod = p1.dot(p2);
	double r1 = p1.norm();
	double r2 = p2.norm();
	double r2_3 = r2 * r2 * r2;
	profile(&eps2, &deps2, r2);
	
	//calculate Legendre function of L and x--cos_angle
	
	int idx1 = int((r1-rBegin) / hStep) - 1;
	int idx2 = int((r2-rBegin) / hStep) - 1; 

	double delta1 = r1 - grid(idx1);
	double delta2 = r2 - grid(idx2);

	
	double u11 = f1(0, idx1) + (f1(0, idx1 + 1) - f1(0, idx1)) * delta1 / hStep;	    								  
	double u12 = f1(0, idx2) + (f1(0, idx2 + 1) - f1(0, idx2)) * delta2 / hStep;	    								    
	double u21 = f2(0, idx1) + (f2(0, idx1 + 1) - f2(0, idx1)) * delta1 / hStep;	    								    
	double u22 = f2(0, idx2) + (f2(0, idx2 + 1) - f2(0, idx2)) * delta2 / hStep;	    								    
	double d11 = f1(1, idx1) + (f1(1, idx1 + 1) - f1(1, idx1)) * delta1 / hStep;	    								    
	double d12 = f1(1, idx2) + (f1(1, idx2 + 1) - f1(1, idx2)) * delta2 / hStep;
	double d21 = f2(1, idx1) + (f2(1, idx1 + 1) - f2(1, idx1)) * delta1 / hStep;	    								    
	double d22 = f2(1, idx2) + (f2(1, idx2 + 1) - f2(1, idx2)) * delta2 / hStep;

	double d2u12 = l * (l + 1) / (r2*r2) - d12 * d12 - (2.0 / r2 + deps2 / eps2) * d12;
	double d2u22 = l * (l + 1) / (r2*r2) - d22 * d22 - (2.0 / r2 + deps2 / eps2) * d22;
	Vector3d g12_deri;
	
	double x1 = p1(0);
	double y1 = p1(1);
	double z1 = p1(2);
	double x2 = p2(0);
	double y2 = p2(1);
	double z2 = p2(2);

	/*
	cout <<"Points " << p1.transpose() << " " << p2.transpose() << endl;
	cout <<"Values " << u11 << " " << u12 << " " << u21 << " " << u22 << endl;
	cout <<"Deriva " << d11 << " " << d12 << " " << d21 << " " << d22 << endl;
	cout <<"DDDDDD " << d2u12 << " " << d2u22 << endl;
	cout <<"Cr12   " << Cr12.transpose() << endl;
	cout <<"g12    " << g12_deri.transpose() << endl;
	*/

	double a = (deps2*r2+2.0*eps2)*plx/(eps2*r2);
	double b = (d12*(d12-d22)+d2u12-d2u22)*plx/(d12-d22); 
	double c = (d22*(d12-d22)+d2u12-d2u22)*plx/(d12-d22); 
	if (r1<r2) {
		double expFact = exp(u11-u12)*(2*l+1);
		g12_deri(0) = -expFact*x2/((d12-d22)*r2_3*eps2) * (a + b - (x1*r2*r2/x2-prod)*dplx/(r1*r2*r2));
		g12_deri(1) = -expFact*y2/((d12-d22)*r2_3*eps2) * (a + b - (y1*r2*r2/y2-prod)*dplx/(r1*r2*r2));
		g12_deri(2) = -expFact*z2/((d12-d22)*r2_3*eps2) * (a + b - (z1*r2*r2/z2-prod)*dplx/(r1*r2*r2));
		if (flagCr12 == 1){
			g12_deri(0) -= (pow(r1/r2,l)*(plx*x2*(-l-1)+dplx*(x1*r2*r2-x2*prod)/(r1*r2))/(r2_3*Cr12(0)));
			g12_deri(1) -= (pow(r1/r2,l)*(plx*y2*(-l-1)+dplx*(y1*r2*r2-y2*prod)/(r1*r2))/(r2_3*Cr12(1)));
			g12_deri(2) -= (pow(r1/r2,l)*(plx*z2*(-l-1)+dplx*(z1*r2*r2-z2*prod)/(r1*r2))/(r2_3*Cr12(2)));
		}
	} else {
		double expFact = exp(u21-u22)*(2*l+1);
		g12_deri(0) = -expFact*x2/((d12-d22)*r2_3*eps2) * (a + c - (x1*r2*r2/x2-prod)*dplx/(r1*r2*r2));
		g12_deri(1) = -expFact*y2/((d12-d22)*r2_3*eps2) * (a + c - (y1*r2*r2/y2-prod)*dplx/(r1*r2*r2));
		g12_deri(2) = -expFact*z2/((d12-d22)*r2_3*eps2) * (a + c - (z1*r2*r2/z2-prod)*dplx/(r1*r2*r2));
		if (flagCr12 == 1){
			g12_deri(0) -= (pow(r2/r1,l)*(plx*x2*l+dplx*(x1*r2*r2-x2*prod)/(r1*r2))/(r1*r2*r2*Cr12(0)));
			g12_deri(1) -= (pow(r2/r1,l)*(plx*y2*l+dplx*(y1*r2*r2-y2*prod)/(r1*r2))/(r1*r2*r2*Cr12(1)));
			g12_deri(2) -= (pow(r2/r1,l)*(plx*z2*l+dplx*(z1*r2*r2-z2*prod)/(r1*r2))/(r1*r2*r2*Cr12(2)));
		} 
	}
	return g12_deri;
}

ostream & DielectricSphere::printObject(ostream &os) {
	os << "Green's Function" << endl;
	os << "Type = Dielectric Sphere" << endl;
	os << "Delta            = " << delta << endl;
	os << "Uniform          = " << uniformFlag << endl;
    os << "Sphere Radius    = " << sphereRadius << endl;    
    os << "Interface Width  = " << interfaceWidth << endl;  
    os << "Eps Inside       = " << epsInside << endl;       
    os << "Eps Outside      = " << epsOutside << endl;      
    os << "Sphere Position  = " << spherePosition.transpose() << endl;  
    os << "Max L Green      = " << maxLGreen << endl;  
    os << "Max L Epsilon    = " << maxLEpsilon << endl;
    os << "h Step           = " << hStep << endl;    
    os << "r Begin          = " << rBegin << endl;   
    os << "r End            = " << rEnd << endl;    
    os << "r Min            = " << rMin << endl;     
    os << "r Max            = " << rMax;
	return os;
}

