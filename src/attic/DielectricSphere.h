#ifndef DIELECTRICSPHERE
#define DIELECTRICSPHERE

typedef Matrix<double, 2, Dynamic> FuncGrid;

class DielectricSphere : public GreensFunction
{
 public:
    DielectricSphere(double epsIn, double epsOut, Vector3d posSph, 
                     double radSph, double widthInt);
    DielectricSphere(Section green);
    ~DielectricSphere(){};
    double evalf(Vector3d &p1, Vector3d &p2);
    double evald(Vector3d &direction, Vector3d &p1, Vector3d &p2);
    double getEpsilon(Vector3d &point);

    /*
    void setHStep(double h){hStep = h;}
    double getHStep(){return hStep;}
    void setRBegin(double r){rBegin = r;}
    double getRBegin(){return rBegin;}
    void setREnd(double r){rEnd = r;}
    double getREnd(){return rEnd;}
    void setRMax(double r){rMax = r;}
    double getRMax(){return rMax;}
    void setRMin(double r){rMin = r;}
    double getRMin(){return rMin;}
    void setMaxLGreen(int l){maxLGreen = l;}
    int getMaxLGreen(){return maxLGreen;}
    void setMaxLEpsilon(int l){maxLEpsilon = l;}
    int getMaxLEpsilon(){return maxLEpsilon;}
    */

    void profile(double *er, double *der_er, double r);
 protected:
    virtual ostream & printObject(std::ostream & os);
 private:    
    void grad(Vector3d &g, Vector3d &p1, Vector3d &p2);
    void initDielectricSphere();
	double getSingleLayer(Vector3d &x, Vector3d &y);
    double getDoubleLayer(Vector3d &p1, Vector3d &p2, Vector3d &n2);
    void independent_solutions();
    void getU1(double r, double r0, FuncGrid & points, int L);
    void getU2(double r, double r0, FuncGrid & points, int L);
    double converged_gf(double r1, double r2, double r12, double cos_theta);
    double greenfunc(double r1, double r2, double r12, double cos_theta, 
                     FuncGrid f1, FuncGrid f2, int l, double epi_r2, 
                     int param_Cr12, double Cr12);
    Vector3d greenfunc_der(Vector3d p1, Vector3d p2, Vector3d Cr12,
                           FuncGrid f1, FuncGrid f2, double plx, 
                           double dplx, int l, int flagCr12);
    Vector3d converged_deri_gf(VectorXd p1, VectorXd p2, double *plx, 
                               double *dplx);
    double sphereRadius;      //central position of the interface
    double interfaceWidth;    //parameter in dertermining teh width of the spherical interface (see jcp 120, 3893 (2004), but should be: D = W/6.0 not W = D/6.0)
    double epsInside;         //permittivity limit inside sphere
    double epsOutside;        //permittivity limit outside sphere
    Vector3d spherePosition;  
    VectorXd grid; //radial solution 1, containing the values and their derivatives
    vector<FuncGrid> radialG1;
    vector<FuncGrid> radialG2;
    FuncGrid radialC1;
    FuncGrid radialC2;
    static int maxLGreen;     //max L in getting the final summation of Green's function
    static int maxLEpsilon;   //max L in getting the constant of C(r1,r2)
    static int nSteps;        //number of integration steps
    static double hStep;      //step size for integration
    static double rBegin;     //start poinr for solving the differential equation in the radial dimension
    static double rEnd;       //final point for solving the differential equation in the radial dimension
    static double rMin;       //min scanning point in radial part in solving differential equation
    static double rMax;       //max scanning point in radial part in solving differential equation

};

#endif

/* Notes

 For now we limit the interfacial variation of epsilon to tanh(x-x0)
 but in theory any profile with reasonable properties should do. 
 The best would be to implement a "profile" class. Not now

 Integration parameters should have defalut values which are class wide

 The dielectric sphere is now located in the origin. SpherePosition is
 therefore not used for now

 */
