#ifndef METALSPHERE
#define METALSPHERE

#include <complex>

typedef complex<double> dcomplex;

class MetalSphere : public GreensFunction<double>
{
 public:
    MetalSphere(double eps, double epsRe, double epsIm, Vector3d &pos,
                double radius);
    MetalSphere(Section green);
    ~MetalSphere(){};
    double evald(Vector3d &direction, Vector3d &p1, Vector3d &p2);
    void compDiagonalElementS(double area);
    void compDiagonalElementD(double area, double radius);
 private:
    double evalGreensFunction(double * source, double * probe);
    double sphRadius;
    double epsSolvent;
    dcomplex epsMetal;
    Vector3d sphPosition;
};
#endif
