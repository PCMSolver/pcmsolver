#ifndef PLANARINTERFACE
#define PLANARINTERFACE

class PlanarInterface : public GreensFunction
{
 public:
    PlanarInterface(double eps1, double eps2, double pos, double width);
    ~PlanarInterface(){};
    double evalf(double* p1, double* p2);
 private:
    double eps1;
    double eps2;
    double pos;
    double width;
    bool computed;
};
#endif
