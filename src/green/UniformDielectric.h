#ifndef UNIFORMDIELECTRIC
#define UNIFORMDIELECTRIC

class UniformDielectric : public GreensFunction
{
 public:
    UniformDielectric(double dielConst);
    UniformDielectric(Section green);
    ~UniformDielectric(){};
    double evalf(Vector3d &p1, Vector3d &p2);
    double evald(Vector3d &direction, Vector3d &p1, Vector3d &p2);
    void setEpsilon(double dielConst) {epsilon = dielConst; };
    double getEpsilon() {return epsilon; };
    double derivative(Vector3d &direction, Vector3d &p1, Vector3d &p2);
    void gradient(Vector3d &gradient, Vector3d &p1, Vector3d &p2);
 private:
    double epsilon;
};
#endif
