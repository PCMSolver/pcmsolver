#ifndef VACUUM
#define VACUUM

template<typename T>
class Vacuum : public GreensFunction<T>
{
 public:
    Vacuum(){GreensFunction<T>::uniformFlag = true;};
    ~Vacuum(){};
    double evald(Eigen::Vector3d & direction, Eigen::Vector3d & p1, Eigen::Vector3d & p2);
    double compDiagonalElementS(double area);
    double compDiagonalElementD(double area, double radius);
    double getDielectricConstant(){return 1.0;}
 private:
    T evalGreensFunction(T * source, T * probe);
};

#endif
