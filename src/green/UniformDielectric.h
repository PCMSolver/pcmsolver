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
    double compDiagonalElementS(double area);
    double compDiagonalElementD(double area, double radius);
 private:
    virtual std::ostream & printObject(std::ostream & os);
    T evalGreensFunction(T * source, T * probe);
    double epsilon;
};
#endif
