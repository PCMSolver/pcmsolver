#ifndef GREENSFUNCTIONSUM
#define GREENSFUNCTIONSUM

class GreensFunctionSum : public GreensFunction
{
 public:
    GreensFunctionSum(GreensFunction &first, GreensFunction &second);
    ~GreensFunctionSum(){
    };
    double evalf(Vector3d &p1, Vector3d &p2);
    double evald(Vector3d &direction, Vector3d &p1, Vector3d &p2, double delta = 0.001);
    double derivative(Vector3d &direction, Vector3d &p1, Vector3d &p2, double delta = 0.001);
    void gradient(Vector3d &gradient, Vector3d &p1, Vector3d &p2, double delta = 0.001);
    GreensFunction* greenFirst;
    GreensFunction* greenSecond;
};
#endif
