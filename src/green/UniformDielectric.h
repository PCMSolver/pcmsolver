#ifndef UNIFORMDIELECTRIC
#define UNIFORMDIELECTRIC

template<class T>
class UniformDielectric : public GreensFunction<T>
{
 public:
    UniformDielectric(double dielConst);
    UniformDielectric(Section green);
    ~UniformDielectric(){};
    double evald(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    void setDielectricConstant(double dielConst) {epsilon = dielConst; };
    double getDielectricConstant(){return epsilon;}
    double compDiagonalElementS(double area);
    double compDiagonalElementD(double area, double radius);
 private:
    virtual std::ostream & printObject(std::ostream & os);
    T evalGreensFunction(T * source, T * probe);
    double epsilon;
};
#endif
