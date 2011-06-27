#ifndef METALSPHERE
#define METALSPHERE

#include <complex>

typedef complex<double> dcomplex;

class MetalSphere : public GreensFunction
{
 public:
    MetalSphere(double eps, double epsRe, double epsIm, Vector3d &pos,
                double radius);
    ~MetalSphere(){};
    double evalf(Vector3d &p1, Vector3d &p2);
    double evald(Vector3d &direction, Vector3d &p1, Vector3d &p2, double delta = 0.001);
 private:
    double sphRadius;
    double epsSolvent;
    dcomplex epsMetal;
    Vector3d sphPosition;
};
#endif
