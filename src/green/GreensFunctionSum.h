#ifndef GREENSFUNCTIONSUM
#define GREENSFUNCTIONSUM

class GreensFunctionSum : public GreensFunction
{
 public:
    GreensFunctionSum(GreensFunction &first, GreensFunction &second);
    GreensFunctionSum(Section green);
    ~GreensFunctionSum(){delete greenFirst; delete greenSecond;}
    double evalf(Vector3d &p1, Vector3d &p2);
    double evald(Vector3d &direction, Vector3d &p1, Vector3d &p2);
    double derivative(Vector3d &direction, Vector3d &p1, Vector3d &p2);
    void gradient(Vector3d &gradient, Vector3d &p1, Vector3d &p2);
 protected:
    GreensFunction* greenFirst;
    GreensFunction* greenSecond;
};
#endif
