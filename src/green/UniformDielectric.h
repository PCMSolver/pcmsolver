#ifndef UNIFORMDIELECTRIC
#define UNIFORMDIELECTRIC

template<class T>
class UniformDielectric : public GreensFunction<T>
{
 public:
    UniformDielectric(double dielConst);
    UniformDielectric(Section green);
    ~UniformDielectric(){};
    double evald(Vector3d &direction, Vector3d &p1, Vector3d &p2);
    void setEpsilon(double dielConst) {epsilon = dielConst; };
    double getEpsilon() {return epsilon; };
 private:
    T evalGreensFunction(T * source, T * probe);
    double epsilon;
};
#endif
